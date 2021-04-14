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
        def proc = ['scp', remote_path.replaceAll(' ', '\\\\ '), dest_path].execute()
        //def proc = ['scp', "\"${remote_path}\"", dest_path].execute()
        def outputStream = new StringBuffer();
        proc.waitForProcessOutput(outputStream, System.err)
        //proc.waitForProcessOutput(System.out, System.err)
        //return outputStream.toString();
        if ((proc.exitValue() != 0) || (!dest.exists())) {
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

def inner_join(ch_a, ch_b) {
    return ch_b.cross(ch_a).map { [it[0][0], *it[1][1..-1], *it[0][1..-1]] }
}
