/*
 * Alignment of one dataset and sequencing type: Metamist cohorts or SGs -> sample sheet -> bwa-mem3 + dupblaster
 * -> CRAM, published with workflow outputs and registered in Metamist. See DESIGN.md.
 */

include { BUILD_SAMPLESHEET } from './modules/local/build_samplesheet'
include { REGISTER_ANALYSIS } from './modules/local/register_analysis'
include { ALIGN } from './subworkflows/local/align'

/*
 * Output paths, relative to outputDir. The output block publishes to these, and REGISTER_ANALYSIS
 * registers the same paths, so the two can't drift apart. Exomes live under exome/, as in cpg-flow.
 */
def typePrefix() {
    params.sequencing_type == 'exome' ? 'exome/' : ''
}

def cramPath(meta) {
    "${typePrefix()}cram/${meta.id}.cram"
}

def metricsPath(meta) {
    "${typePrefix()}qc/markduplicates_metrics/${meta.id}.markduplicates-metrics"
}

// per-run records: the sample sheet, the CRAM index, and any SGs that failed
def runPath(name) {
    "alignment_samplesheets/${params.arGuid}/${name}"
}

def metamistProject() {
    params.accessLevel == 'test' ? "${params.dataset}-test" : params.dataset
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
}

/*
 * Outputs this session has passed the pre-flight check for. -resume keeps the session ID and needs the same
 * work dir, so a resumed run finds the claims of its earlier attempts here.
 */
def claimsFile() {
    workflow.workDir.resolve("nf-align-claims/${workflow.sessionId}.txt")
}

/*
 * Fail before any alignment if an output this run would publish already exists. overwrite = 'false'
 * would skip such a file silently, and registration would then record a file this run didn't write.
 * Metamist-built sheets leave out SGs with a registered CRAM, so a clash is most likely an unregistered CRAM.
 * On -resume, outputs claimed by an earlier attempt of this session are expected, and are republished.
 */
def checkNoExistingOutputs(rows) {
    def outputs = rows
        .collect { meta, _r -> meta }
        .unique { meta -> meta.id }
        .collectMany { meta -> [cramPath(meta), "${cramPath(meta)}.crai", metricsPath(meta)] }
        .collect { rel -> uriOf(workflow.outputDir.resolve(rel)) }
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
    main:
    validateParams()
    def kinds = ['fastq', 'fastq_ora', 'cram', 'bam', 'aligned']
    // matches {genome,exome}_realignment_shards in cpg-flow's config_template.toml
    // CLI params arrive as strings, and `1..'2'` would silently range over char codes
    def shards = params.realignment_shards.toInteger() ?: (params.sequencing_type == 'genome' ? 6 : 1)

    // the FASTA and everything that sits beside it: its .fai, and the bwa-mem3 index built with the FASTA as prefix.
    // Staged side by side, so the FASTA's own name is the bwa-mem3 index prefix. A .dict is deliberately left out (RISKS.md)
    def ref_suffixes = ['fai', 'bwt.2bit.64', 'amb', 'ann', 'pac']
    def reference = channel.value([file(params.ref_fasta), ref_suffixes.collect { ext -> file("${params.ref_fasta}.${ext}") }])

    // there are three entrypoints: --cohorts, --sgids and --samplesheet. If a samplesheet is given as input, it
    // shouldn't be written to the output folder. built_sheet vs. sheet differentiates the two cases
    def built_sheet = channel.empty()
    def sheet = channel.empty()
    if (params.samplesheet) {
        sheet = channel.fromPath(params.samplesheet, checkIfExists: true)
    }
    else {
        // one of the two is empty: SGs are discovered through cohorts, or queried directly
        def cohorts = params.cohorts ? idList(params.cohorts) : []
        def sgids = params.sgids ? idList(params.sgids) : []
        built_sheet = BUILD_SAMPLESHEET(metamistProject(), params.sequencing_type, cohorts, sgids)
        sheet = built_sheet
    }

    def rows = sheet
        .splitCsv(header: true, sep: '\t')
        .map { r ->
            if (!(r.kind in kinds)) {
                fail("Unknown kind '${r.kind}' for ${r.sg_id}, expected one of ${kinds}")
            }
            def meta = [id: r.sg_id, kind: r.kind, units: r.n_units.toInteger(), ref: r.reads_ref ?: params.ref_fasta]
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
    if (params.register) {
        registered = REGISTER_ANALYSIS(
            ALIGN.out.cram.map { meta, _cram, _crai -> [meta, uriOf(workflow.outputDir.resolve(cramPath(meta)))] },
            metamistProject(),
            params.sequencing_type,
            params.arGuid,
            params.register_wait_minutes.toInteger(),
        )
    }

    // with errorStrategy 'ignore' a failed SG doesn't fail the run, so name each one and the step it failed at
    def failed = checked.map { meta, _r -> [meta.id, 'expected'] }
        .mix(ALIGN.out.cram.map { meta, _cram, _crai -> [meta.id, 'aligned'] })
        .mix(registered.map { meta -> [meta.id, 'registered'] })
        .groupTuple()
        .map { id, events ->
            def failed_at = !('aligned' in events) ? 'align' : (params.register && !('registered' in events) ? 'register' : null)
            [id: id, failed_at: failed_at]
        }
        .filter { record -> record.failed_at }
    failed.toList().subscribe { records ->
        if (records) {
            log.warn("${records.size()} SG(s) failed and were skipped (see ${runPath('failed.csv')}): " +
                records.collect { record -> "${record.id} (${record.failed_at})" }.join(', '))
        }
    }

    publish:
    samplesheet = built_sheet
    crams = ALIGN.out.cram.map { meta, cram, crai -> [id: meta.id, cram: cram, crai: crai] }
    markdup = ALIGN.out.markdup.map { meta, metrics -> [id: meta.id, metrics: metrics] }
    failed = failed
}

// paths are relative to outputDir, which defaults to the dataset bucket (nextflow.config)
output {
    samplesheet {
        path { sheet -> sheet >> runPath('samplesheet.tsv') }
    }
    crams {
        path { record ->
            record.cram >> cramPath(record)
            record.crai >> "${cramPath(record)}.crai"
        }
        index {
            path runPath('crams.csv')
            header true
        }
    }
    markdup {
        path { record -> record.metrics >> metricsPath(record) }
    }
    failed {
        index {
            path runPath('failed.csv')
            header true
        }
    }
}
