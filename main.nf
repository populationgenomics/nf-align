/*
 * Alignment of one sequencing type: Metamist cohorts or SGs, from any number of datasets -> sample sheet ->
 * bwa-mem3 + dupblaster -> CRAM, published to and registered in each SG's own dataset. Run records go to the
 * analysis dataset (params.dataset). See DESIGN.md.
 */

include { BUILD_SAMPLESHEET } from './modules/local/build_samplesheet'
include { REGISTER_ANALYSIS } from './modules/local/register_analysis'
include { WRITE_RUN_RECORDS } from './modules/local/write_run_records'
include { ALIGN } from './subworkflows/local/align'

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

// each file ALIGN_TO_CRAM and MERGE_TO_CRAM write, by name -> where it goes under the SG's output root.
// Plain String keys: a GString key would never match the filename publishDir looks up
def publishPaths(id) {
    ["cram/${id}.cram", "cram/${id}.cram.crai", "qc/markduplicates_metrics/${id}.markduplicates-metrics"]
        .collectEntries { rel -> [(rel.tokenize('/').last()): "${typePrefix()}${rel}".toString()] }
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
    def missing = ['dataset', 'sequencing_type', 'arGuid', 'ref_fasta'].findAll { name -> !params[name] }
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
    if (!(params.align_container && params.cpg_container)) {
        error('The workflow needs --align_container and --cpg_container')
    }
    if (!(params.register.toString() in ['true', 'false'])) {
        error("--register must be true or false, got '${params.register}'")
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
 * Fail before any alignment if an output this run would publish already exists, in any SG's dataset.
 * publishDir's overwrite: false leaves such a file in place silently, and registration would then record a
 * file this run didn't write (ADR 0009).
 * Metamist-built sheets leave out SGs with a registered CRAM, so a clash is most likely an unregistered CRAM.
 * On -resume, outputs claimed by an earlier attempt of this session are expected, and are republished.
 */
def checkNoExistingOutputs(rows) {
    def outputs = rows
        .collect { meta, _r -> meta }
        .unique { meta -> meta.id }
        .collectMany { meta -> meta.publish.keySet().collect { filename -> destination(meta, filename) } }
    def claims = claimsFile()
    def claimed = claims.exists() ? claims.readLines().toSet() : [] as Set
    def existing = outputs
        .findAll { uri -> !(uri in claimed) }
        .findAll { uri -> file(uri).exists() }
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

    // the FASTA and everything that sits beside it: its .fai, and the bwa-mem3 index built with the FASTA as prefix.
    // Staged side by side, so the FASTA's own name is the bwa-mem3 index prefix. A .dict is deliberately left out (RISKS.md)
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
            def meta = [
                id: r.sg_id,
                kind: r.kind,
                units: r.n_units.toInteger(),
                ref: r.reads_ref ?: params.ref_fasta,
                dataset: r.dataset,
                project: metamistProject(r.dataset),
                outRoot: datasetRoot(r.dataset),
                publish: publishPaths(r.sg_id),
            ]
            [meta, r]
        }
        .branch { meta, _r ->
            aligned: meta.kind == 'aligned'
            to_align: true
        }

    // hold every row until the whole sheet has been checked, so a clash fails the run before any work starts
    // this scenario covers alignment outputs which exist in GCS but haven't been registered. We don't want to attempt
    // regenaration/overwriting, but they would appear in the samplesheet based on the metamist-only logic
    def checked = rows.to_align
        .toList()
        .flatMap { items -> checkNoExistingOutputs(items) }

    ALIGN(checked, rows.aligned, reference, shards, params.dupblaster_extra_args)

    def registered = channel.empty()
    if (register) {
        registered = REGISTER_ANALYSIS(
            ALIGN.out.cram.map { meta, _cram, _crai -> [meta, destination(meta, "${meta.id}.cram")] },
            params.sequencing_type,
            params.arGuid,
            params.ref_fasta,
            params.register_wait_minutes.toInteger(),
        )
    }

    // with errorStrategy 'ignore' a failed SG doesn't fail the run, so name each one and the step it failed at
    def failed = checked.map { meta, _r -> [meta.id, meta.dataset, 'expected'] }
        .mix(ALIGN.out.cram.map { meta, _cram, _crai -> [meta.id, meta.dataset, 'aligned'] })
        .mix(registered.map { meta -> [meta.id, meta.dataset, 'registered'] })
        .groupTuple(by: [0, 1])
        .map { id, dataset, events ->
            def failed_at = !('aligned' in events) ? 'align' : (register && !('registered' in events) ? 'register' : null)
            [sg_id: id, dataset: dataset, failed_at: failed_at]
        }
        .filter { record -> record.failed_at }
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
