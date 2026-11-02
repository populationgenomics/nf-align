/*
 * Alignment of one sequencing type: Metamist cohorts or SGs, from any number of datasets -> sample sheet ->
 * bwa-mem3 + dupblaster -> CRAM, fingerprinted and QC'd, published to and registered in each SG's own dataset. Run
 * records go to the analysis dataset (params.dataset). See DESIGN.md.
 */

include { BUILD_SAMPLESHEET } from './modules/local/build_samplesheet'
include { REGISTER_ANALYSIS as REGISTER_CRAM } from './modules/local/register_analysis'
include { REGISTER_ANALYSIS as REGISTER_CRAM_QC } from './modules/local/register_analysis'
include { WRITE_RUN_RECORDS } from './modules/local/write_run_records'
include { SOMALIER_EXTRACT } from './modules/local/somalier_extract'
include { ALIGN } from './subworkflows/local/align'
include { CRAM_QC } from './subworkflows/local/cram_qc'

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

// per-run records, in the analysis dataset: the sample sheet, the CRAM index, and any SGs that failed
def runDir() {
    "${datasetRoot(params.dataset)}/alignment_samplesheets/${params.arGuid}"
}

// `--cram_qc false` arrives as the string 'false', which Groovy treats as true
def cramQcOn() {
    params.cram_qc.toString().toBoolean()
}

// an SG's QC metrics, by their secondary-file key in its cram_qc analysis -> file name (ADR 0014). The duplicate
// metrics come from alignment, so they're published even when CRAM QC is off
def cramQcFiles(id) {
    def coverage = params.sequencing_type == 'genome' ? 'wgs' : 'hs'
    def qc = !cramQcOn() ? [:] : [
        verifybamid: "${id}.verify-bamid.selfSM",
        samtools_stats: "${id}.samtools-stats",
        alignment_summary_metrics: "${id}.alignment_summary_metrics",
        base_distribution_by_cycle_metrics: "${id}.base_distribution_by_cycle_metrics",
        insert_size_metrics: "${id}.insert_size_metrics",
        quality_by_cycle_metrics: "${id}.quality_by_cycle_metrics",
        quality_yield_metrics: "${id}.quality_yield_metrics",
        ("picard_${coverage}_metrics".toString()): "${id}.picard-${coverage}-metrics",
    ]
    (qc + [markduplicates_metrics: "${id}.markduplicates-metrics"]).collectEntries { key, name -> [(key): name.toString()] }
}

// each file the SG's tasks publish, by name -> where it goes under the SG's output root: the CRAM and its
// companions in cram/, and every QC file in one cram_qc/<sg>/ folder (ADR 0014).
// Plain String keys: a GString key would never match the filename publishDir looks up
def publishPaths(id) {
    def cram = ["${id}.cram", "${id}.cram.crai", "${id}.cram.somalier"].collect { name -> "cram/${name}" }
    def qc = cramQcFiles(id).values().collect { name -> "cram_qc/${id}/${name}" }
    (cram + qc).collectEntries { rel -> [(rel.tokenize('/').last()): "${typePrefix()}${rel}".toString()] }
}

// the one place a published URI is built from meta
def destination(meta, filename) {
    "${meta.outRoot}/${meta.publish[filename.toString()]}".toString()
}

// gs://bucket/path for cloud paths (Path.toString() drops the scheme), a plain path locally
def uriOf(path) {
    path.scheme == 'file' ? path.toString() : path.toUriString()
}

// a list from config, or a comma/space separated string from the command line
def idList(param) {
    (param instanceof List) ? param : param.toString().tokenize(', ')
}

def validateParams() {
    def missing = ['dataset', 'sequencing_type', 'arGuid', 'ref_fasta', 'aligned_ref_fasta', 'somalier_sites'].findAll { name -> !params[name] }
    if (missing) {
        error("Missing required params: ${missing.collect { name -> '--' + name }.join(', ')}")
    }
    if (!(params.sequencing_type in ['genome', 'exome'])) {
        error("--sequencing_type must be genome or exome, got '${params.sequencing_type}'")
    }
    if (!(params.accessLevel in ['full', 'standard', 'test'])) {
        error("--accessLevel must be full, standard or test, got '${params.accessLevel}'")
    }
    // a bare `--dupblaster_extra_args`, or `--dupblaster_extra_args ''`, arrives as 'true'
    if (!(params.dupblaster_extra_args instanceof CharSequence) || params.dupblaster_extra_args.toString() == 'true') {
        error("--dupblaster_extra_args must be a string of dupblaster arguments, got '${params.dupblaster_extra_args}'")
    }
    if ([params.cohorts, params.sgids, params.samplesheet].count { input -> input } != 1) {
        error('Provide exactly one of --cohorts, --sgids or --samplesheet')
    }
    if (!(params.align_container && params.cpg_container && params.somalier_container)) {
        error('The workflow needs --align_container, --cpg_container and --somalier_container')
    }
    if (!(params.register.toString() in ['true', 'false'])) {
        error("--register must be true or false, got '${params.register}'")
    }
    if (!(params.cram_qc.toString() in ['true', 'false'])) {
        error("--cram_qc must be true or false, got '${params.cram_qc}'")
    }
    if (cramQcOn()) {
        if (!(params.verifybamid_container && params.picard_container)) {
            error('CRAM QC needs --verifybamid_container and --picard_container (or --cram_qc false)')
        }
        def qc_missing = cramQcRefParams().values().findAll { name -> !params[name] }
        if (qc_missing) {
            error("CRAM QC needs ${qc_missing.collect { name -> '--' + name }.join(', ')} (or --cram_qc false)")
        }
    }
}

// the param names of the references CRAM QC reads, for the run's sequencing type: the alignment reference's
// .dict, VerifyBamID's contamination resources, and Picard's interval list
def cramQcRefParams() {
    def contamination = ['ud', 'bed', 'mu'].collectEntries { ext -> [(ext): "${params.sequencing_type}_contam_${ext}".toString()] }
    def intervals = params.sequencing_type == 'genome' ? 'genome_coverage_interval_list' : 'exome_evaluation_interval_list'
    [dict: 'ref_fasta_dict'] + contamination + [intervals: intervals]
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
    def existing = ['samplesheet.tsv', 'crams.csv', 'failed.csv']
        .collect { name -> "${runDir()}/${name}" }
        .findAll { uri -> file(uri).exists() }
    if (existing) {
        error(
            "Run records for arGuid ${params.arGuid} already exist. Launch with a new --arGuid:\n  ${existing.join('\n  ')}"
        )
    }
}

/*
 * Outputs this session has passed the pre-flight check for. -resume keeps the session ID and needs the same
 * work dir, so a resumed run finds the claims of its earlier attempts here.
 */
def claimsFile() {
    workflow.workDir.resolve("nf-align-claims/${workflow.sessionId}.txt")
}

/*
 * Fail before any alignment or fingerprinting if an output this run would publish already exists, in any SG's
 * dataset. Already-aligned SGs publish nothing, so they never get here.
 * publishDir's overwrite: false leaves such a file in place silently, and registration would then record a
 * file this run didn't write. The run only ever registers what it wrote (ADR 0009, ADR 0012).
 * Metamist-built sheets leave out SGs with a registered CRAM, so a clash is most likely an unregistered CRAM.
 * On -resume, outputs claimed by an earlier attempt of this session are expected, and are republished.
 */
def checkNoExistingOutputs(rows) {
    def claims = claimsFile()
    def claimed = claims.exists() ? claims.readLines().toSet() : [] as Set

    def outputs = rows
        .collect { meta, _r -> meta }
        .unique { meta -> meta.id }
        .collectMany { meta -> meta.publish.keySet().collect { filename -> destination(meta, filename) } }
    def existing = outputs.findAll { uri -> existsUnclaimed(uri, claimed) }
    if (existing) {
        fail(
            "${existing.size()} output(s) already exist, and are probably not registered in Metamist. " +
            "Register or delete them, then rerun:\n  ${existing.join('\n  ')}"
        )
    }
    claims.parent.mkdirs()
    claims.text = (claimed + outputs).join('\n')
    return rows
}

// an existing output that no earlier attempt of this session published
def existsUnclaimed(uri, claimed) {
    !(uri in claimed) && file(uri).exists()
}

// error() inside a channel operator reaches the console only as "Unexpected error", so log the reason first
def fail(message) {
    log.error(message)
    error(message)
}

workflow {
    validateParams()
    checkNoRunRecords()
    def kinds = ['fastq', 'fastq_ora', 'cram', 'bam', 'aligned']
    // matches {genome,exome}_realignment_shards in cpg-flow's config_template.toml
    // CLI params arrive as strings, and `1..'2'` would silently range over char codes
    def shards = params.realignment_shards.toInteger() ?: (params.sequencing_type == 'genome' ? 6 : 1)
    // likewise `--register false` arrives as the string 'false', which Groovy treats as true
    def register = params.register.toString().toBoolean()
    def cram_qc = cramQcOn()

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
        sheet = BUILD_SAMPLESHEET(params.accessLevel, params.sequencing_type, cohorts, sgids, runDir())
    }

    def rows = sheet
        .splitCsv(header: true, sep: '\t')
        .map { r ->
            if (!(r.kind in kinds)) {
                fail("Unknown kind '${r.kind}' for ${r.sg_id}, expected one of ${kinds}")
            }
            // a CRAM decodes only against the reference it was written with. A BAM's ref is staged but never read
            if (r.kind == 'cram' && !r.reads_ref) {
                fail("${r.sg_id}: a cram row needs reads_ref, the reference the CRAM was written against")
            }
            // decides where the SG is published and registered, and which dataset its compute is billed to
            if (!r.dataset) {
                fail("${r.sg_id}: a row needs dataset, the dataset the SG belongs to")
            }
            def aligned = r.kind == 'aligned'
            def meta = [
                id: r.sg_id,
                kind: r.kind,
                units: r.n_units.toInteger(),
                // an aligned CRAM whose analysis doesn't record its reference predates nf-align (ADR 0006)
                ref: r.reads_ref ?: (aligned ? params.aligned_ref_fasta : params.ref_fasta),
                dataset: r.dataset,
                project: metamistProject(r.dataset),
                outRoot: datasetRoot(r.dataset),
                // an aligned SG's CRAM is already registered, fingerprint or not, so nothing is published
                publish: aligned ? [:] : publishPaths(r.sg_id),
            ]
            [meta, r]
        }
        .branch { meta, _r ->
            aligned: meta.kind == 'aligned'
            to_align: true
        }

    // hold every row until the whole sheet has been checked, so a clash fails the run before any work starts
    // this scenario covers alignment outputs which exist in GCS but haven't been registered. We don't want to attempt
    // regenaration/overwriting, or to register files this run didn't write, so they fail the run
    def checked = rows.to_align
        .toList()
        .flatMap { items -> checkNoExistingOutputs(items) }

    ALIGN(checked, rows.aligned, reference, shards, params.dupblaster_extra_args)

    // fingerprint each newly aligned CRAM against ref_fasta, which it was written with. Only the FASTA and .fai
    // are staged, not the bwa-mem3 index
    def fingerprinted = SOMALIER_EXTRACT(
        ALIGN.out.cram.map { meta, cram, crai -> [meta + [ref: params.ref_fasta], cram, crai, file(params.ref_fasta), file("${params.ref_fasta}.fai")] },
        file(params.somalier_sites),
    )

    // QC each newly aligned CRAM, alongside its fingerprint and apart from its registration: a QC failure never
    // holds up the CRAM (ADR 0014). Only the FASTA, .fai and .dict are staged
    def qc = channel.empty()
    if (cram_qc) {
        def refs = cramQcRefParams().collectEntries { key, name -> [(key): file(params[name])] }
        qc = CRAM_QC(
            ALIGN.out.cram,
            channel.value([file(params.ref_fasta), file("${params.ref_fasta}.fai"), refs.dict]),
            channel.value([refs.ud, refs.bed, refs.mu]),
            channel.value(refs.intervals),
            params.sequencing_type,
        )
    }

    // one analysis per SG: the CRAM, with its index and fingerprint as secondary files. A CRAM whose fingerprint
    // failed isn't registered (ADR 0012)
    def wait_minutes = params.register_wait_minutes.toInteger()
    def registered = channel.empty()
    def qc_registered = channel.empty()
    if (register) {
        registered = REGISTER_CRAM(
            fingerprinted.map { meta, _somalier ->
                def companions = [index: "${meta.id}.cram.crai", somalier: "${meta.id}.cram.somalier"]
                def analysis_meta = [sequencing_type: params.sequencing_type, 'ar-guid': params.arGuid, ref_fasta: meta.ref]
                [meta, destination(meta, "${meta.id}.cram"), companions.collectEntries { key, name -> [(key): destination(meta, name)] }, analysis_meta]
            },
            'cram',
            wait_minutes,
        )

        // and one cram_qc analysis per SG: the CRAM again, with every QC file as a secondary file. Only once the
        // CRAM's own analysis exists, so it never points at an unregistered CRAM (ADR 0014)
        qc_registered = REGISTER_CRAM_QC(
            qc.map { meta, _files -> [meta.id, meta] }
                .join(registered.map { meta -> [meta.id, true] })
                .map { _id, meta, _registered ->
                    def files = cramQcFiles(meta.id).collectEntries { key, name -> [(key): destination(meta, name)] }
                    [meta, destination(meta, "${meta.id}.cram"), files, [sequencing_type: params.sequencing_type, 'ar-guid': params.arGuid]]
                },
            'cram_qc',
            wait_minutes,
        )
    }

    // with errorStrategy 'ignore' a failed SG doesn't fail the run, so name each one and the first step it didn't
    // get through: each SG is aligned, fingerprinted, registered, QC'd and has its QC registered
    // I'm not convinced by this whole mechanic, it's a lot of complexity for something we're able to check in the UI
    def steps = [aligned: 'align', fingerprinted: 'somalier'] +
        (register ? [registered: 'register'] : [:]) +
        (cram_qc ? [qc_done: 'qc'] : [:]) +
        (register && cram_qc ? [qc_registered: 'register_qc'] : [:])
    def failed = checked.map { meta, _r -> [meta.id, meta.dataset, 'to_align'] }
        .mix(ALIGN.out.cram.map { meta, _cram, _crai -> [meta.id, meta.dataset, 'aligned'] })
        .mix(fingerprinted.map { meta, _somalier -> [meta.id, meta.dataset, 'fingerprinted'] })
        .mix(registered.map { meta -> [meta.id, meta.dataset, 'registered'] })
        .mix(qc.map { meta, _files -> [meta.id, meta.dataset, 'qc_done'] })
        .mix(qc_registered.map { meta -> [meta.id, meta.dataset, 'qc_registered'] })
        .groupTuple(by: [0, 1])
        .flatMap { id, dataset, events ->
            def failed_at = steps.find { event, _step -> !(event in events) }?.value
            failed_at ? [[sg_id: id, dataset: dataset, failed_at: failed_at]] : []
        }
    failed.toList().subscribe { records ->
        if (records) {
            log.warn("${records.size()} SG(s) failed and were skipped (see ${runDir()}/failed.csv): " +
                records.collect { record -> "${record.sg_id} (${record.failed_at})" }.join(', '))
        }
    }

    // the published URIs of each newly aligned CRAM, as worked out in meta
    def crams = ALIGN.out.cram.map { meta, _cram, _crai ->
        [sg_id: meta.id, dataset: meta.dataset, cram: destination(meta, "${meta.id}.cram"), crai: destination(meta, "${meta.id}.cram.crai")]
    }
    // sorted, so the records don't depend on task completion order and -resume finds the task cached.
    // Each emits [] when nothing arrives, so both files are always written
    def by_sg = { a, b -> a.sg_id <=> b.sg_id }
    WRITE_RUN_RECORDS(crams.toSortedList(by_sg), failed.toSortedList(by_sg), runDir())
}
