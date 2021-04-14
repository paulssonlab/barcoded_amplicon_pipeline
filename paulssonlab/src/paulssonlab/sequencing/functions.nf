import nextflow.util.CsvParser

def scp(remote_path, dest_path) {
    def dest = file(dest_path)
    if (!dest.exists()) {
        def dest_parent = dest.getParent()
        if (!dest_parent.exists()) {
            dest_parent.mkdirs()
        }
        // SEE: https://stackoverflow.com/questions/25300550/difference-in-collecting-output-of-executing-external-command-in-groovy
        // and https://stackoverflow.com/questions/159148/groovy-executing-shell-commands
        // without escaping spaces, would need -T to deal with scp quoted source paths
        // SEE: https://stackoverflow.com/questions/54598689/scp-fails-with-protocol-error-filename-does-not-match-request
        // not sure why I need two backslashes to escape spaces
        def command = ['scp', remote_path.replaceAll(' ', '\\\\ '), dest_path]
        def proc = command.execute()
        //def proc = ['scp', "\"${remote_path}\"", dest_path].execute()
        def outputStream = new StringBuffer();
        proc.waitForProcessOutput(outputStream, System.err)
        //proc.waitForProcessOutput(System.out, System.err)
        //return outputStream.toString()
        if ((proc.exitValue() != 0) || (!dest.exists())) {
            println "scp failed with output:"
            println outputStream.toString()
            throw new Exception("scp of '${remote_path}' to '${dest_path}' failed")
        }
    }
    return dest
}

def read_tsv(path) {
    def parser = new CsvParser()
        .setSeparator('\t')
    def data = []
    def tsv_file = file(path)
    tsv_file.eachLine { line ->
        data << parser.parse(line)
    }
    return data
}

def join_key(ch_a, ch_b, old_key, new_key) {
    return ch_b.cross(ch_a.map { [it[old_key], it] } ).map { [(new_key): it[0][1], *:it[1][1]] }
}

def join_each(ch_a, ch_b, old_key, new_key) {
    return ch_a.combine(ch_b).map { a, b ->
        [(new_key): a[old_key].collect { b.get(it) }, *:a]
    }
}
