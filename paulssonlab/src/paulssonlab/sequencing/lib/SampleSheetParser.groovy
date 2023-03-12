@Grab(group='com.moandjiezana.toml', module='toml4j', version='0.7.2')
import com.moandjiezana.toml.Toml

import nextflow.util.CsvParser
import java.nio.file.Paths
import java.math.BigDecimal

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

    private static List load(String path, Map defaults = [:], substitute = true) {
        def sampleSheet = new Toml().read(new File(path)).toMap()
        def samples = []
        def tsv = sampleSheet.get("tsv")
        if (tsv) {
            def tsvParser = new CsvParser()
                .setSeparator('\t')
            def data = []
            def header
            tsv.eachLine { line ->
                def parsedLine = tsvParser.parse(line)
                parsedLine = parsedLine.collect { it =~ /\s*-\s*/ ? "" : it }
                if (!header) {
                    header = parsedLine
                } else {
                    if (parsedLine.size() != header.size()) {
                        throw new Exception("Expecting ${header.size()} columns in tsv, got ${parsedLine.size()}")
                    }
                    samples << [header, parsedLine].transpose().collectEntries()
                }
            }
        }
        def samplesTable = sampleSheet.get("samples")
        samplesTable?.each {
            samples << it
        }
        def paramSets = sampleSheet.getOrDefault("params", [[run_path: "default"]])
        paramSets = paramSets.collect {
            [run_path: "default", *:it]
        }
        if (anyDuplicates(paramSets*.run_path)) {
            throw new Exception("Param sets must have unique run paths")
        }
        def runs = paramSets.collectMany { p ->
            samples.collect { s -> [*:defaults, *:p, *:s] }
        }
        def engine = new groovy.text.SimpleTemplateEngine()
        runs.eachWithIndex { it, index ->
            it.id = index
            def meta = it.clone() // so substitutions can't depend on each other
            it.replaceAll { k, v ->
                try {
                    v = new BigDecimal(v)
                } catch (NumberFormatException e) {}
                if (v instanceof String && (substitute == true || (substitute instanceof Collection && k in substitute))) {
                    engine.createTemplate(v).make(meta).toString()
                } else {
                    v
                }
            }
            it.run_path = Paths.get(it.run_path, it["name"]) as String
        }
        if (anyDuplicates(runs*.run_path)) {
            throw new Exception("Runs must have unique run paths")
        }
        return runs
    }
}
