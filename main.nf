/*
 * Alignment of one sequencing type: Metamist cohorts or SGs, from any number of datasets -> sample sheet ->
 * bwa-mem3 + dupblaster -> CRAM, fingerprinted and QC'd, published to and registered in each SG's own dataset. Run
 * records go to the analysis dataset (params.dataset). With --qc_aligned, already-aligned SGs without CRAM QC are
 * QC'd too, without realigning (ADR 0015). See DESIGN.md.
 */

include { BUILD_SAMPLESHEET } from './modules/local/build_samplesheet'
include { PREFLIGHT } from './modules/local/preflight'
include { REGISTER_ANALYSIS as REGISTER_CRAM } from './modules/local/register_analysis'
include { REGISTER_ANALYSIS as REGISTER_CRAM_QC } from './modules/local/register_analysis'
include { WRITE_RUN_RECORDS } from './modules/local/write_run_records'
include { SOMALIER_EXTRACT } from './modules/local/somalier_extract'
include { ALIGN } from './subworkflows/local/align'
include { CRAM_QC } from './subworkflows/local/cram_qc'

/*
 * Types only: the defaults are in nextflow.config, which conf/base.config reads before this script is parsed.
 * Nextflow converts command-line values to these types (`--cram_qc false` is a Boolean), and fails the run if a
 * param without `?` has no value.
 */
params {
    dataset: String
    sequencing_type: String
    accessLevel: String
    arGuid: String
    cohorts: String? // or a list, from a config file
    sgids: String?   // or a list, from a config file
    samplesheet: Path?
    output_root: String?

    register: Boolean
    cram_qc: Boolean
    qc_aligned: Boolean
    register_wait_minutes: Integer
    realignment_shards: Integer
    dupblaster_extra_args: String

    ref_fasta: String
    ref_fasta_dict: String
    aligned_ref_fasta: String
    somalier_sites: String

    align_container: String
    cpg_container: String
    somalier_container: String

    // CRAM QC's, needed only with --cram_qc, and the references only for the run's sequencing type
    verifybamid_container: String?
    picard_container: String?
    genome_contam_ud: String?
    genome_contam_bed: String?
    genome_contam_mu: String?
    exome_contam_ud: String?
    exome_contam_bed: String?
    exome_contam_mu: String?
    genome_coverage_interval_list: String?
    exome_evaluation_interval_list: String?
}

/*
 * Output destinations. Each SG's are worked out here and carried in its meta: publishDir (conf/base.config),
 * the pre-flight check and registration all read them from there, so they can't drift apart.
 * Exomes live under exome/, as in cpg-flow.
 */
def typePrefix() {
    params.sequencing_type == 'exome' ? 'exome/' : ''
}

def metamistProject(ds) {
    params.accessLevel == 'test' ? "${ds}-test" : ds
}

// a dataset's bucket, or <output_root>/<ds> for local runs (absolute, whatever the launch dir)
def datasetRoot(ds) {
    params.output_root
        ? uriOf(file(params.output_root).resolve(ds))
        : "gs://cpg-${ds}-${params.accessLevel == 'test' ? 'test' : 'main'}"
}

// per-run records, in the analysis dataset: the sample sheet and the CRAM index
def runDir() {
    "${datasetRoot(params.dataset)}/alignment_samplesheets/${params.arGuid}"
}

// the files CRAM QC's own tools write for an SG, by their secondary-file key in its cram_qc analysis -> file name
// (ADR 0014)
def qcToolFiles(id) {
    def coverage = params.sequencing_type == 'genome' ? 'wgs' : 'hs'
    def suffixes = [
        verifybamid: 'verify-bamid.selfSM',
        samtools_stats: 'samtools-stats',
        alignment_summary_metrics: 'alignment_summary_metrics',
        base_distribution_by_cycle_metrics: 'base_distribution_by_cycle_metrics',
        insert_size_metrics: 'insert_size_metrics',
        quality_by_cycle_metrics: 'quality_by_cycle_metrics',
        quality_yield_metrics: 'quality_yield_metrics',
        ("picard_${coverage}_metrics".toString()): "picard-${coverage}-metrics",
    ]
    suffixes.collectEntries { key, suffix -> [(key): "${id}.${suffix}".toString()] }
}

// the files of an SG's cram analysis (ADR 0012): the CRAM, then its secondary files
def cramFiles(meta) {
    [cram: meta.cram, index: "${meta.cram}.crai".toString(), somalier: "${meta.cram}.somalier".toString()]
}

// the QC files this run publishes for an SG, the secondary files of its cram_qc analysis (ADR 0014): its QC tools'
// files (with --cram_qc), then the duplicate metrics, which come from alignment, so a QC-only SG has none (ADR 0015)
def qcFiles(meta) {
    def names = (params.cram_qc ? qcToolFiles(meta.id) : [:]) +
        (meta.qcOnly ? [:] : [markduplicates_metrics: "${meta.id}.markduplicates-metrics".toString()])
    names.collectEntries { key, name -> [(key): "${meta.qcDir}/${name}".toString()] }
}

// every file this run publishes for an SG, which the pre-flight check makes sure don't exist yet. An aligned SG's
// CRAM is already registered, so at most its QC is published
def outputs(meta) {
    if (meta.kind != 'aligned') {
        return cramFiles(meta) + qcFiles(meta)
    }
    meta.qcOnly ? qcFiles(meta) : [:]
}

// a QC-only SG's duplicate metrics, if the run that aligned it published them in the SG's dataset: nf-align's
// cram_qc/<sg>/, or cpg-flow's qc/markduplicates_metrics/. [markduplicates_metrics: uri], or [:] if neither
// exists. Not this run's file, so it's only looked for, never published or checked (ADR 0015)
def existingDuplicateMetrics(meta) {
    def name = "${meta.id}.markduplicates-metrics"
    def found = ["${meta.qcDir}/${name}", "${datasetRoot(meta.dataset)}/${typePrefix()}qc/markduplicates_metrics/${name}"]
        .collect { uri -> uri.toString() }
        .find { uri -> file(uri).exists() }
    found ? [markduplicates_metrics: found] : [:]
}

// the .dict of the reference a CRAM was written against: ref_fasta's is its own param, and any other FASTA's
// (an already-aligned CRAM's) is the <name>.dict beside <name>.fasta, where Picard looks for it. null for a
// FASTA without a FASTA extension
def refDict(fasta) {
    if (fasta == params.ref_fasta) {
        return params.ref_fasta_dict
    }
    def dict = fasta.replaceFirst(/\.(fasta|fa|fna)(\.gz)?$/, '.dict')
    dict == fasta ? null : dict
}

// gs://bucket/path for cloud paths (Path.toString() drops the scheme), a plain path locally
def uriOf(path) {
    path.scheme == 'file' ? path.toString() : path.toUriString()
}

// a list from config, or a comma/space separated string from the command line
def idList(param) {
    (param instanceof List) ? param : param.toString().tokenize(', ')
}

// the param names of the references CRAM QC reads, for the run's sequencing type: VerifyBamID's contamination
// resources, and Picard's interval list
def cramQcRefParams() {
    def type = params.sequencing_type
    ['ud', 'bed', 'mu'].collectEntries { ext -> [(ext): "${type}_contam_${ext}".toString()] } +
        [intervals: type == 'genome' ? 'genome_coverage_interval_list' : 'exome_evaluation_interval_list']
}

// what the params block can't check: allowed values, and params that are needed together or not at all
def validateParams() {
    if (!(params.sequencing_type in ['genome', 'exome'])) {
        error("--sequencing_type must be genome or exome, got '${params.sequencing_type}'")
    }
    if (!(params.accessLevel in ['full', 'standard', 'test'])) {
        error("--accessLevel must be full, standard or test, got '${params.accessLevel}'")
    }
    // a bare `--dupblaster_extra_args`, or `--dupblaster_extra_args ''`, arrives as 'true'
    if (params.dupblaster_extra_args == 'true') {
        error("--dupblaster_extra_args must be a string of dupblaster arguments, got '${params.dupblaster_extra_args}'")
    }
    if ([params.cohorts, params.sgids, params.samplesheet].count { input -> input } != 1) {
        error('Provide exactly one of --cohorts, --sgids or --samplesheet')
    }
    if (params.qc_aligned && !params.cram_qc) {
        error('--qc_aligned needs CRAM QC, so it can\'t be combined with --cram_qc false')
    }
    def qc_params = ['verifybamid_container', 'picard_container'] + cramQcRefParams().values()
    def qc_missing = params.cram_qc ? qc_params.findAll { name -> !params[name] } : []
    if (qc_missing) {
        error("CRAM QC needs ${qc_missing.collect { name -> '--' + name }.join(', ')} (or --cram_qc false)")
    }
}

/*
 * Fail before the sample sheet is built if this arGuid's run records already exist: publishDir would keep the
 * old ones silently. arGuid is unique in production, but set by hand while testing. A resumed run finds its
 * own records there, so -resume skips the check (ADR 0009).
 */
def checkNoRunRecords() {
    if (workflow.resume) {
        return
    }
    def existing = ['samplesheet.tsv', 'crams.csv']
        .collect { name -> "${runDir()}/${name}" }
        .findAll { uri -> file(uri).exists() }
    if (existing) {
        error(
            "Run records for arGuid ${params.arGuid} already exist. Launch with a new --arGuid:\n  ${existing.join('\n  ')}"
        )
    }
}

/*
 * A sample sheet row -> [meta, row]. meta holds the SG's destinations, and ref, the reference its CRAM is (or will
 * be) written against. The row keeps what only alignment reads, such as the input reads' own reference.
 */
def parseRow(r) {
    def kinds = ['fastq', 'fastq_ora', 'cram', 'bam', 'aligned']
    if (!(r.kind in kinds)) {
        fail("Unknown kind '${r.kind}' for ${r.sg_id}, expected one of ${kinds}")
    }
    // a CRAM decodes only against the reference it was written with. The sample sheet builder fills it in for
    // an aligned SG whose analysis doesn't record one (ADR 0006)
    if (r.kind in ['cram', 'aligned'] && !r.reads_ref) {
        fail("${r.sg_id}: ${r.kind} rows need reads_ref, the reference the CRAM was written against")
    }
    // decides where the SG is published and registered, and which dataset its compute is billed to
    if (!r.dataset) {
        fail("${r.sg_id}: a row needs dataset, the dataset the SG belongs to")
    }
    def aligned = r.kind == 'aligned'
    def base = "${datasetRoot(r.dataset)}/${typePrefix()}"
    def meta = [
        id: r.sg_id,
        kind: r.kind,
        units: r.n_units.toInteger(),
        dataset: r.dataset,
        project: metamistProject(r.dataset),
        ref: aligned ? r.reads_ref : params.ref_fasta,
        // an aligned SG's registered CRAM, or where this run publishes the SG's new one
        cram: aligned ? r.reads : "${base}cram/${r.sg_id}.cram".toString(),
        cramDir: "${base}cram".toString(),
        qcDir: "${base}cram_qc/${r.sg_id}".toString(),
        // with --qc_aligned, an already-aligned SG whose CRAM has no cram_qc analysis gets CRAM QC, and nothing
        // else (ADR 0015). A --samplesheet without the cram_qc column has every aligned SG QC'd
        qcOnly: aligned && params.qc_aligned && !r.cram_qc,
    ]
    if (meta.qcOnly && !refDict(meta.ref)) {
        fail("${r.sg_id}: can't QC a CRAM written against ${meta.ref}, as Picard needs the .dict beside a .fasta, .fa or .fna")
    }
    [meta, r]
}

// error() inside a channel operator reaches the console only as "Unexpected error", so log the reason first
def fail(message) {
    log.error(message)
    error(message)
}

workflow {
    validateParams()
    checkNoRunRecords()
    // matches {genome,exome}_realignment_shards in cpg-flow's config_template.toml
    def shards = params.realignment_shards ?: (params.sequencing_type == 'genome' ? 6 : 1)

    // the FASTA and everything that sits beside it: its .fai, the bwa-mem3 index built with the FASTA as prefix,
    // and the .alt that makes alignment alt-aware. Staged side by side, so the FASTA's own name is the bwa-mem3
    // index prefix. The alt-annotated .dict goes too: bwa-mem3 copies its @SQ lines, AH tags and all (ADR 0010)
    def ref_suffixes = ['fai', 'bwt.2bit.64', 'amb', 'ann', 'pac', 'alt']
    def reference = channel.value([file(params.ref_fasta), file(params.ref_fasta_dict), ref_suffixes.collect { ext -> file("${params.ref_fasta}.${ext}") }])

    // there are three entrypoints: --cohorts, --sgids and --samplesheet. Only a sheet built here is published
    // to the run records; a --samplesheet input isn't
    def sheet = channel.empty()
    if (params.samplesheet) {
        sheet = channel.fromPath(params.samplesheet, checkIfExists: true)
    }
    else {
        // one of the two is empty: SGs are discovered through cohorts, or queried directly
        def cohorts = params.cohorts ? idList(params.cohorts) : []
        def sgids = params.sgids ? idList(params.sgids) : []
        sheet = BUILD_SAMPLESHEET(params.accessLevel, params.sequencing_type, cohorts, sgids, params.aligned_ref_fasta, runDir())
    }

    // hold every row until the whole sheet has passed the pre-flight check, so a clash fails the run before any
    // work starts. Already-aligned SGs that aren't QC-only publish nothing, and go no further
    def rows = PREFLIGHT(
        sheet
            .splitCsv(header: true, sep: '\t')
            .map { r -> parseRow(r) }
            .toList()
            .map { items -> [items, items.collectMany { meta, _r -> outputs(meta).values() }.unique()] }
    )
        .flatMap { items -> items }
        .branch { meta, _r ->
            to_align: meta.kind != 'aligned'
            qc_only: meta.qcOnly
        }

    ALIGN(rows.to_align, reference, shards, params.dupblaster_extra_args)

    // fingerprint each newly aligned CRAM against ref_fasta, which it was written with. Only the FASTA and .fai
    // are staged, not the bwa-mem3 index
    def fingerprinted = SOMALIER_EXTRACT(
        ALIGN.out.cram.map { meta, cram, crai -> [meta, cram, crai, file(meta.ref), file("${meta.ref}.fai")] },
        file(params.somalier_sites),
    )

    // QC each newly aligned CRAM, alongside its fingerprint and apart from its registration: a QC failure never
    // holds up the CRAM (ADR 0014). And each QC-only SG's registered CRAM (ADR 0015). Each against the reference it
    // was written with: only the FASTA, .fai and .dict are staged
    def qc = channel.empty()
    if (params.cram_qc) {
        def refs = cramQcRefParams().collectEntries { key, name -> [(key): file(params[name])] }
        qc = CRAM_QC(
            ALIGN.out.cram
                .mix(rows.qc_only.map { meta, r -> [meta, file(r.reads), file(r.reads_index)] })
                .map { meta, cram, crai -> [meta, cram, crai, file(meta.ref), file("${meta.ref}.fai"), file(refDict(meta.ref))] },
            channel.value([refs.ud, refs.bed, refs.mu]),
            channel.value(refs.intervals),
            params.sequencing_type,
        )
    }

    if (params.register) {
        // one analysis per SG: the CRAM, with its index and fingerprint as secondary files. A CRAM whose fingerprint
        // failed isn't registered (ADR 0012)
        def registered = REGISTER_CRAM(
            fingerprinted.map { meta, _somalier ->
                def files = cramFiles(meta)
                [meta, files.cram, files.subMap('index', 'somalier'), [sequencing_type: params.sequencing_type, 'ar-guid': params.arGuid, ref_fasta: meta.ref]]
            },
            'cram',
            params.register_wait_minutes,
        )

        // and one cram_qc analysis per SG: the CRAM again, with every QC file as a secondary file. Only once the
        // CRAM's own analysis exists, so it never points at an unregistered CRAM (ADR 0014). A QC-only SG's is
        // registered already, and its duplicate metrics are whatever its aligning run left, if anything (ADR 0015)
        def has_cram_analysis = registered.mix(rows.qc_only.map { meta, _r -> meta })
        REGISTER_CRAM_QC(
            qc.map { meta, _files -> [meta.id, meta] }
                .join(has_cram_analysis.map { meta -> [meta.id, true] })
                .map { _id, meta, _registered ->
                    def files = qcFiles(meta) + (meta.qcOnly ? existingDuplicateMetrics(meta) : [:])
                    [meta, meta.cram, files, [sequencing_type: params.sequencing_type, 'ar-guid': params.arGuid]]
                },
            'cram_qc',
            params.register_wait_minutes,
        )
    }

    // the published URIs of each newly aligned CRAM, sorted so the records don't depend on task completion order and
    // -resume finds the task cached. Emits [] when nothing arrives, so the file is always written
    def crams = ALIGN.out.cram.map { meta, _cram, _crai ->
        [sg_id: meta.id, dataset: meta.dataset, cram: meta.cram, crai: "${meta.cram}.crai".toString()]
    }
    WRITE_RUN_RECORDS(crams.toSortedList { a, b -> a.sg_id <=> b.sg_id }, runDir())
}
