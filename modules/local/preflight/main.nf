/*
 * Fail the run before any work starts if a file it would publish already exists. publishDir's
 * overwrite: false would leave such a file in place silently, and registration would then record a file this run
 * didn't write.
 * A task, not a function in main.nf, so that -resume skips it: once it has passed, Nextflow caches it, and a resumed
 * run doesn't trip over the files its own earlier attempt published. One that failed isn't cached, so it checks
 * again.
 * Runs natively in the Nextflow process (exec:), so it needs no container or VM.
 */
process PREFLIGHT {
    input:
    tuple val(rows), val(outputs) // the sample sheet rows, held back until checked; every URI the run would publish

    output:
    val rows

    exec:
    def existing = outputs.findAll { uri -> file(uri).exists() }
    if (existing) {
        def message = "${existing.size()} output(s) already exist, and are probably not registered in Metamist. " +
            "Register or delete them, then rerun:\n  ${existing.join('\n  ')}"
        log.error(message)
        error(message)
    }
}
