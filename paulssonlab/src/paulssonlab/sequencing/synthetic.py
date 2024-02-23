import itertools as it
import tempfile

import numpy as np


def mutagenize_seq(seq, q=0, error=0, letters="ATCG", rng=None):
    # mark errors as upper-case to make debugging easier (GraphAligner doesn't care)
    letters = list(letters)
    if rng is None:
        rng = np.random.default_rng()
    if q and error:
        raise ValueError("at most one of q and error can be specified")
    if q:
        error = 10 ** (-q / 10)
    num_errors = rng.binomial(len(seq), error)
    error_indices = rng.choice(len(seq), size=num_errors)
    for idx in error_indices:
        seq = seq[:idx] + rng.choice(letters) + seq[idx + 1 :]
    return seq


def generate_reads(segments, num_reads=100, q=0, error=0, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    num_choices = np.array([len(s) for s in segments])
    num_segments = len(segments)
    true_path = rng.integers(num_choices[np.newaxis, :], size=(num_reads, num_segments))
    reversed = rng.integers(2, size=num_reads)
    reads = []
    for read_idx in range(num_reads):
        read = "".join(
            [
                variants[variant_idx]
                for variants, variant_idx in zip(segments, true_path[read_idx])
            ]
        )
        read = mutagenize_seq(read, q=q, error=error, rng=rng)
        ###### TODO
        # read = (
        #     read[0] + mutagenize_seq(read[1:-1], q=q, error=error, rng=rng) + read[-1]
        # )
        ######
        if reversed[read_idx]:
            read = str(sequence.reverse_complement(read))
        ####### TODO
        # read = "N" + read
        # read = read[:-1] + "G"
        # read = "G" + read[1:]
        #######
        reads.append(read)
    # add trailing newline
    formatted_reads = (
        "\n".join([f">r{idx}\n{read}" for idx, read in enumerate(reads)]) + "\n"
    )
    ground_truth = dict(true_path=true_path, reversed=reversed)
    return formatted_reads, ground_truth


def generate_gfa(segments):
    lines = ["H\tVN:Z:1.0"]
    lines.extend(
        [
            f"S\ts{s}={v}\t{seq}"
            for s, variants in enumerate(segments)
            for v, seq in enumerate(variants)
        ]
    )
    lines.extend(
        [
            f"L\ts{s}={v1}\t+\ts{s+1}={v2}\t+\t0M"
            for s in range(len(segments) - 1)
            for v1, v2 in it.product(
                range(len(segments[s])), range(len(segments[s + 1]))
            )
        ]
    )
    return "\n".join(lines) + "\n"  # add trailing newline


def run_aligner(gfa_filename, reads_filename, args=[]):
    cmd_base = ["/home/jqs1/micromamba/envs/graphaligner/bin/GraphAligner"]
    # cmd_base = ["/home/jqs1/paulsson-home/bin/GraphAligner"]
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".gaf") as gaf_file:
        cmd = [
            *cmd_base,
            "-g",
            gfa_filename,
            "-f",
            reads_filename,
            "-a",
            gaf_file.name,
            *args,
        ]
        start = time.time()
        out = subprocess.run(cmd, capture_output=True)
        stop = time.time()
        if out.returncode != 0:
            print("STDOUT:")
            print(out.stdout.decode())
            print()
            print("STDERR:")
            print(out.stderr.decode())
            print()
            raise RuntimeError("GraphAligner returned non-zero exit status")
        runtime = stop - start
        print("STDOUT")
        print(out.stdout.decode())
        print("STDERR")
        print(out.stderr.decode())
        gaf = sio.read_gaf(gaf_file.name)
        return gaf, runtime


def run_aligner_synthetic(segments, args=[["-x", "vg"]], num_reads=4, q=0, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    with (
        tempfile.NamedTemporaryFile(mode="w+", suffix=".gfa") as gfa_file,
        tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta") as reads_file,
    ):
        gfa = generate_gfa(segments)
        # print(gfa)
        gfa_file.write(gfa)
        gfa_file.flush()
        reads, ground_truth = generate_reads(
            segments, num_reads=num_reads, q=q, rng=rng
        )
        # print(reads)
        reads_file.write(reads)
        reads_file.flush()
        res = []
        for cmd_args in args:
            res.append(run_aligner(gfa_file.name, reads_file.name, args=cmd_args))
    return res, ground_truth


def check_path_equality(path, true_path):
    if path[0][0] == "<":
        path = path[::-1]
    if len(path) != len(true_path):
        return False
    for segment_idx, p in enumerate(path):
        match = re.match(r"(?:<|>)s(\d+)=(\d+)", p)
        if int(match.group(1)) != segment_idx:
            return False
        if int(match.group(2)) != true_path[segment_idx]:
            return False
    return True


def check_alignment(gaf, ground_truth):
    errors = set()
    for idx in range(len(gaf)):
        path = gaf.column("path")[idx].as_py()
        if not check_path_equality(path, ground_truth["true_path"][idx]):
            # TODO
            # print(">>>",path,ground_truth["true_path"][idx])
            errors.add(idx)
        if (path[0][0] == "<") != ground_truth["reversed"][idx]:
            errors.add(idx)
    return errors
