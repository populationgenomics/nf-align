// crams.csv for the run records directory, from the records main.nf collected.
// Runs natively in the Nextflow process (exec:), so it needs no container or VM.
process WRITE_RUN_RECORDS {
    input:
    val crams   // [sg_id, dataset, cram, crai] maps, one per newly aligned SG
    val run_dir // only read by publishDir (conf/base.config)

    output:
    path 'crams.csv'

    exec:
    task.workDir.resolve('crams.csv').text = csvText(['sg_id', 'dataset', 'cram', 'crai'], crams)
}

// values are SG IDs, dataset names and URIs, which never contain commas, so no quoting is needed
def csvText(columns, records) {
    ([columns] + records.collect { record -> columns.collect { column -> record[column] } })
        .collect { row -> row.join(',') + '\n' }
        .join('')
}
