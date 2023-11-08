import re
from enum import IntEnum
from itertools import takewhile

from gfapy import Gfa

from paulssonlab.util.sequence import reverse_complement

# FROM: https://github.com/jeffdaily/parasail/blob/600fb26151ff19899ee39a214972dcf2b9b11ed7/src/cigar.c#L18
# and
# FROM: https://github.com/kcleal/pywfa
# for enum trickery, see https://www.notinventedhere.org/articles/python/how-to-use-strings-as-name-aliases-in-python-enums.html
CigarOp = IntEnum(
    "CigarOp",
    [
        ("M", 0),
        ("I", 1),
        ("D", 2),
        ("N", 3),
        ("S", 4),
        ("H", 5),
        ("P", 6),
        ("=", 7),
        ("X", 8),
        ("B", 9),
    ],
)
CigarOp.__repr__ = lambda self: self.name
CigarOp.__str__ = CigarOp.__repr__
CigarOp.__format__ = lambda self, spec: self.__repr__()

OP_TO_COLUMN_NAME = {
    CigarOp["="]: "matches",
    CigarOp.X: "mismatches",
    CigarOp.I: "insertions",
    CigarOp.D: "deletions",
}


def encode_cigar(cigar):
    return "".join(f"{length}{op}" for op, length in cigar)


def decode_cigar(s):
    return [
        (CigarOp[match[2]], int(match[1]))
        for match in re.finditer(r"(\d+)(M|I|D|N|S|H|P|=|X|B)", s)
    ]


def _parse_variant(s):
    try:
        return int(s)
    except:
        return s


def cut_cigar(
    cigar,
    path,
    name_to_seq,
    sequence=None,
    variant_sep="=",
    key_sep="|",
    return_sequences=True,
    return_indices=True,
    return_counts=True,
    return_cigars=True,
    cigar_as_string=True,
    return_variants=True,
    separate_ends=True,
):
    if isinstance(name_to_seq, Gfa):
        name_to_seq = sgfa.gfa_name_mapping(name_to_seq)
    segment_lengths = [len(name_to_seq[name]) for name in path]
    segment_names = [name[1:] for name in path]
    segment_rc = [name[0] == "<" for name in path]
    # if separate_ends:
    #     segment_names = ["upstream", *segment_names, "downstream"]
    #     segment_lengths = [0, *segment_lengths, 0]
    #     segment_rc = [False, *segment_rc, False]
    ops = [c[0] for c in cigar]
    op_lengths = [c[1] for c in cigar]
    segment_idx = 0
    cigar_idx = 0
    query_idx = 0
    segment_length = segment_lengths[segment_idx]
    segment_name = segment_names[segment_idx]
    if "=" in segment_name:
        eq_idx = segment_name.index("=")
        segment_name, variant_name = segment_name[:eq_idx], _parse_variant(
            segment_name[eq_idx + 1 :]
        )
    else:
        variant_name = None
    res = {}
    if separate_ends:
        res["upstream"] = {}
        res["downstream"] = {}
        upstream_insertions = list(takewhile(lambda x: x[0] == CigarOp.I, cigar))
        downstream_insertions = list(
            reversed(list(takewhile(lambda x: x[0] == CigarOp.I, reversed(cigar))))
        )
        cigar = cigar[
            len(upstream_insertions) : len(cigar) - len(downstream_insertions)
        ]
        upstream_length = sum(x[1] for x in upstream_insertions)
        downstream_length = sum(x[1] for x in downstream_insertions)
        if return_sequences:
            if sequence is not None:
                res["upstream"]["seq"] = sequence[:upstream_length]
                res["downstream"]["seq"] = sequence[len(sequence) - downstream_length :]
        if return_indices:
            pass
        if return_counts:
            res["upstream"] = {}
        if return_cigars:
            pass
    res[segment_name] = {}
    if variant_name is not None:
        res[segment_name]["variant"] = variant_name
    if return_indices:
        res[segment_name]["start"] = query_idx
    if return_counts:
        for col_name in OP_TO_COLUMN_NAME.values():
            res[segment_name][col_name] = 0
    if return_cigars:
        res[segment_name]["cigar"] = []
    op = ops[cigar_idx]
    op_length = op_lengths[cigar_idx]
    while True:
        advance = min(x for x in (op_length, segment_length) if x is not None)
        print("!", segment_idx, advance)
        # print(
        #     f"op {op} {op_length} seg {segment_name} {segment_length} advance {advance}"
        # )
        if op in [CigarOp.I, CigarOp["="], CigarOp.X]:
            query_idx += advance
        op_length -= advance
        if op in [CigarOp.D, CigarOp["="], CigarOp.X]:
            segment_length -= advance
        if return_counts:
            col_name = OP_TO_COLUMN_NAME[op]
            res[segment_name][col_name] += advance
        if return_cigars:
            res[segment_name]["cigar"].append((op, advance))
        if return_variants:
            pass
        print("*", segment_length, op_length, f"{segment_idx}/{len(segment_lengths)}")
        if segment_length == 0:
            if return_indices:
                res[segment_name]["end"] = query_idx
                res[segment_name]["reverse_complement"] = segment_rc[segment_idx]
            if return_sequences:
                if sequence is not None:
                    segment_seq = sequence[
                        res[segment_name]["start"] : res[segment_name]["end"]
                    ]
                    if segment_rc[segment_idx]:
                        segment_seq = reverse_complement(segment_seq)
                    res[segment_name]["seq"] = segment_seq
            if return_cigars:
                if segment_rc[segment_idx]:
                    res[segment_name]["cigar"].reverse()
                if cigar_as_string:
                    print(">", segment_name)
                    res[segment_name]["cigar"] = encode_cigar(
                        res[segment_name]["cigar"]
                    )
            segment_idx += 1
            if segment_idx == len(segment_lengths):
                segment_length = None
            else:
                segment_length = segment_lengths[segment_idx]
                segment_name = segment_names[segment_idx]
                if "=" in segment_name:
                    eq_idx = segment_name.index("=")
                    segment_name, variant_name = (
                        segment_name[:eq_idx],
                        _parse_variant(segment_name[eq_idx + 1 :]),
                    )
                else:
                    variant_name = None
                res[segment_name] = {}
                if variant_name is not None:
                    res[segment_name]["variant"] = variant_name
                if return_indices:
                    res[segment_name]["start"] = query_idx
                if return_counts:
                    for col_name in OP_TO_COLUMN_NAME.values():
                        res[segment_name][col_name] = 0
                if return_cigars:
                    res[segment_name]["cigar"] = []
        if op_length == 0:
            cigar_idx += 1
            if cigar_idx == len(ops):
                print("HELP")
                pass  # can we ever get here without immediately breaking below?
            else:
                op = ops[cigar_idx]
                op_length = op_lengths[cigar_idx]
        # TODO: need to wrap up
        if cigar_idx == len(ops) and segment_idx == len(segment_lengths):
            break
    if key_sep is not None:
        row = {}
        for segment_name, segment_info in res.items():
            for k, v in segment_info.items():
                row[f"{segment_name}{key_sep}{k}"] = v
        res = row
    return res
