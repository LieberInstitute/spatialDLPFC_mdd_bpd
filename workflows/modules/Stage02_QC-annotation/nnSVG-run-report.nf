// Run on the launcher after publishing finishes, including when sample tasks fail.
// No R installation is needed on the launcher, and no failed-task output channel
// is required. Work-directory diagnostics also remain available for inspection.
def writeNnsvgRunReport(report_path, trace_path, check_slurm = false) {
    def report = file(report_path)
    if (!report.exists()) {
        return // Preparation did not finish; there is no sample manifest to summarize.
    }
    def columns = ['sample_id', 'nnSVG_seed', 'n_rows', 'min.genes_per.spot',
                   'min.spots_per.gene', 'exit_status', 'error_desc']
    def rows = report.readLines().drop(1).findAll { it }.collect { line ->
        // Preparation writes only validated sample IDs, integer seeds, and placeholders.
        def fields = line.split(',', -1).collect { it.replaceAll('^"|"$', '') }
        [sample_id: fields[0], nnSVG_seed: fields[1], n_rows: 'NA',
         'min.genes_per.spot': 'NA', 'min.spots_per.gene': 'NA',
         exit_status: 'NA', error_desc: 'Not run']
    }
    def tasks = [:]
    def preparation_complete = false
    def trace = file(trace_path)
    if (trace.exists()) {
        def lines = trace.readLines()
        def header = lines[0].split('\t', -1).toList()
        lines.drop(1).findAll { it }.each { line ->
            def fields = line.split('\t', -1)
            def task = [:]
            header.eachWithIndex { key, i -> task[key] = fields[i] }
            if (task.process.tokenize(':').last() == 'NNSVG_PREP' && task.status in ['COMPLETED', 'CACHED']) {
                preparation_complete = true
            }
            if (task.process.tokenize(':').last() == 'NNSVG_PER_SAMPLE') {
                tasks[task.tag] = task
            }
        }
    }
    if (!preparation_complete) {
        return // Do not alter a previous run's report when this run never prepared samples.
    }
    rows.each { row ->
        def task = tasks[row.sample_id]
        if (task != null) {
            def work_dir = file(task.workdir)
            def metrics = work_dir.resolve('nnSVG_run_specifics.tsv')
            if (metrics.exists()) {
                def lines = metrics.readLines()
                if (lines.size() > 1) {
                    def header = lines[0].split('\t', -1)
                    def fields = lines[1].split('\t', -1)
                    header.eachWithIndex { key, i -> row[key] = fields[i] }
                }
            }
            // The task exit code is authoritative, including container/scheduler failures.
            row.exit_status = task.exit in ['-', '2147483647'] ? 'NA' : task.exit
            if (task.status in ['COMPLETED', 'CACHED'] && task.exit == '0') {
                row.error_desc = ''
            } else {
                def log_lines = ['.command.err', '.command.log'].collectMany { name ->
                    def log_file = work_dir.resolve(name)
                    log_file.exists() ? log_file.readLines() : []
                }
                def timeout_lines = log_lines.findAll { line ->
                    line.toLowerCase() =~ /due to time limit|timed out|timeout|time limit exceeded|elapsed time limit|reached elapsed time limit/
                }.unique()
                // Slurm accounting can identify TIMEOUT/OOM even when R was killed
                // before writing an error. Query only failed HPC tasks, with a bound.
                def scheduler_state = ''
                if (check_slurm && task.native_id && task.native_id != '-') {
                    try {
                        def accounting = new ProcessBuilder(['sacct', '-n', '-X', '-j', task.native_id,
                            '--format=State,ExitCode', '--parsable2']).redirectErrorStream(true).start()
                        if (accounting.waitFor(10, java.util.concurrent.TimeUnit.SECONDS)) {
                            if (accounting.exitValue() == 0) {
                                scheduler_state = accounting.inputStream.text.trim()
                            }
                        } else {
                            accounting.destroyForcibly()
                        }
                    } catch (Exception ignored) {
                        // Accounting may be unavailable; retain the task logs and exit code.
                    }
                }
                if (scheduler_state =~ /(?i)TIMEOUT/) {
                    row.error_desc = 'Timeout: Slurm reports ' + scheduler_state.replaceAll(/[\t\r\n]+/, ' ')
                } else if (timeout_lines) {
                    row.error_desc = 'Timeout: ' + timeout_lines.join(' ').replaceAll(/[\t\r\n]+/, ' ')
                } else if (scheduler_state =~ /(?i)OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|CANCELLED/) {
                    row.error_desc = 'Slurm reports ' + scheduler_state.replaceAll(/[\t\r\n]+/, ' ')
                } else if (row.error_desc in ['', 'Not completed', 'Not run']) {
                    def error_lines = log_lines.findAll { line ->
                        line =~ /(?i)error|killed|out.of.memory|oom|cancelled|failed/
                    }.unique()
                    row.error_desc = error_lines ? error_lines.join(' ').replaceAll(/[\t\r\n]+/, ' ') :
                        "Task ${task.status} (exit ${task.exit}); no R error recorded. See ${task.workdir}"
                }
            }
        }
    }
    def csv_lines = [columns.join(',')]
    rows.each { row ->
        csv_lines.add(columns.collect { key ->
            '"' + row[key].toString().replace('"', '""') + '"'
        }.join(','))
    }
    // Replace the published manifest, never the cached preparation input.
    def temporary = report.resolveSibling(report.fileName.toString() + '.tmp')
    temporary.text = csv_lines.join('\n') + '\n'
    java.nio.file.Files.move(temporary, report, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
}
