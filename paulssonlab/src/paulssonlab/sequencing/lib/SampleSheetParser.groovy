@Grab(group='com.moandjiezana.toml', module='toml4j', version='0.7.2')
import com.moandjiezana.toml.Toml

import nextflow.util.CsvParser
import java.nio.file.Paths

class SampleSheetParser {
    private static Boolean anyDuplicates(x) {
        return x.toUnique().size() != x.size()
    }

    private static def renameKey(map, oldKey, newKey, defaultValue) {
        return map.put(newKey, map.remove(oldKey) ?: defaultValue)
    }

    private static def splitString(str) {
        if (str == null) {
            return []
        }
        else if (str instanceof CharSequence) {
            return str.split(",").collect { it.strip() }
        }
        else {
            return str
        }
    }

    private static List load(String path) {
        def sampleSheet = new Toml().read(new File(path)).toMap()
        def samples = []
        def tsv = sampleSheet.get("tsv")
        if (tsv) {
            def tsvParser = new CsvParser()
                .setSeparator('\t')
            def data = []
            tsv.eachLine { line ->
                def parsedLine = tsvParser.parse(line)
                if (parsedLine.size() != 2) {
                    throw new Exception("Expecting two columns in tsv, got ${parsedLine.size()}")
                }
                samples << [["reads_prefix", "references"], parsedLine].transpose().collectEntries()
            }
        }
        def samplesTable = sampleSheet.get("samples")
        if (tsv && samplesTable) {
            throw new Exception("Cannot specify both tsv and samples in sample sheet TOML")
        }
        samplesTable?.each {
            samples << it
        }
        samples.each {
            it.reference_names = splitString(it.get("references"))
            it.remove("references")
            if (!it.get("name")) {
                it.name = it.getOrDefault("reads_prefix", "default")
            }
        }
        if (anyDuplicates(samples*.getOrDefault("name", ""))) {
            throw new Exception("Samples must have unique names")
        }
        def paramSets = sampleSheet.getOrDefault("params", [[name: "default"]])
        paramSets.each { renameKey(it, "name", "param_set", "default") }
        if (anyDuplicates(paramSets*.param_set)) {
            throw new Exception("Param sets must have unique names")
        }
        def runs = paramSets.collectMany { p ->
            samples.collect { s -> [*:p, *:s] }
        }
        runs.eachWithIndex { it, index ->
            it.id = index
            it.run_path = Paths.get(it.param_set, it.get("name")) as String
        }
        return runs
    }
}
