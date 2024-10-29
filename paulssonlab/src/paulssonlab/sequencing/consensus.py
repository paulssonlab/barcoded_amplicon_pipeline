import hashlib

import awkward as ak
import numba
import numpy as np
import pyarrow as pa

from paulssonlab.util.sequence import reverse_complement

try:
    import pyabpoa
except ImportError:
    pass

try:
    import spoa
except ImportError:
    pass

GAP_CHAR = ord("-")
SPACE_CHAR = ord(" ")

# SEE: https://github.com/yangao07/abPOA/tree/main/python
ABPOA_DEFAULTS = {"aln_mode": "l"}
# SEE: https://github.com/rvaser/spoa
# AND https://github.com/nanoporetech/pyspoa/blob/master/pyspoa.cpp
SPOA_DEFAULTS = {"algorithm": 0}


# FROM: https://github.com/yangao07/abPOA/blob/c6c9dacf92414cd1a358cbd1d28b6f2ca37b30ed/src/abpoa_output.c#L270
def abpoa_coverage_to_phred(coverage, num_seqs):
    coverage = np.asarray(coverage)
    x = 13.8 * (1.25 * coverage / num_seqs - 0.25)
    p = 1 - 1.0 / (1.0 + np.exp(-x))
    return (-10 * np.log10(p) + 0.499).astype(np.uint8)


def poa(
    seqs,
    phreds=None,
    method="abpoa",
    num_consensus_seqs=1,
    return_msa=False,
    return_phreds=False,
    min_frequency=0.25,
    sort=True,
    **kwargs,
):
    if sort:
        # sort reads by length, see https://github.com/yangao07/abPOA/issues/48
        # this is what medaka does for abpoa
        sort_idxs = sorted(
            range(len(seqs)), key=lambda idx: len(seqs[idx]), reverse=True
        )
        seqs = [seqs[idx] for idx in sort_idxs]
        if phreds is not None:
            phreds = list(phreds[sort_idxs])
    if method == "abpoa":
        aligner = pyabpoa.msa_aligner(**{**ABPOA_DEFAULTS, **kwargs})
        res = aligner.msa(
            seqs,
            # TODO: need to implement phred input
            # weights=phreds,
            out_cons=num_consensus_seqs > 0,
            out_msa=return_msa,
            max_n_cons=num_consensus_seqs or 0,
            min_freq=min_frequency,
        )
        consensus_seqs = res.cons_seq
        num_seqs = len(seqs)
        if return_phreds:
            consensus_phreds = [
                abpoa_coverage_to_phred(cov, num_seqs) for cov in res.cons_cov
            ]
        msa_seqs = res.msa_seq
    elif method == "spoa":
        if num_consensus_seqs > 1:
            raise ValueError("spoa can only output 1 consensus sequence")
        if phreds is not None:
            # TODO: would be easy to implement
            raise ValueError(
                "spoa does support weighting by phred scores but that is not currently implemented by pyspoa"
            )
        if return_phreds:
            raise ValueError("spoa does not support outputting consensus phred scores")
        consensus_seq, msa_seqs = spoa.poa(
            list(seqs), genmsa=return_msa, **{**SPOA_DEFAULTS, **kwargs}
        )
        consensus_seqs = [consensus_seq]
    elif method == "first":
        # return first sequence, used only for debugging
        consensus_seqs = [seqs[0]]
        consensus_phreds = [phreds[0]] if phreds else None
        msa_seqs = None
    else:
        raise ValueError(f"invalid method {method}")
    res = dict(consensus_seqs=consensus_seqs)
    if return_phreds:
        res["consensus_phreds"] = consensus_phreds
    if return_msa:
        res["msa_seqs"] = msa_seqs
    return res


def prepare_reads(seqs, rcs, phreds=None):
    if isinstance(seqs, pa.Array):
        seqs = seqs.to_pylist()  # prefer seqs to be given as a python list
    if isinstance(rcs, pa.Array):
        rcs = ak.from_arrow(rcs)
    seqs_oriented = [
        reverse_complement(seq) if rc else seq for seq, rc in zip(seqs, rcs)
    ]
    res = dict(seqs=seqs_oriented)
    if phreds is not None:
        phreds_oriented = ak.from_arrow(
            pa.array(
                [
                    phred.values.to_numpy()[::-1] if rc else phred.values.to_numpy()
                    for phred, rc in zip(phreds, rcs)
                ]
            )
        )
        res["phreds"] = phreds_oriented
    return res


def _names_to_hash(names, algorithm="sha256"):
    h = hashlib.new(algorithm)
    for name in names:
        h.update(name.encode())
    return h.hexdigest()


def _get_consensus(
    seqs,
    rc,
    phreds=None,
    names=None,
    return_phreds=True,
    **kwargs,
):
    prepared_reads = prepare_reads(
        seqs,
        rc,
        phreds=phreds,
    )
    if len(seqs) == 0:
        raise ValueError("getting the consensus of an empty list of sequences")
    elif len(seqs) == 1:
        res = dict(consensus_seq=seqs[0])
        if phreds is not None and return_phreds:
            res["consensus_phred"] = phreds[0]
        if names is not None:
            res["name"] = names[0]
        return res
    poa_result = poa(
        prepared_reads["seqs"],
        phreds=prepared_reads.get("phreds"),
        num_consensus_seqs=1,
        return_msa=False,
        return_phreds=return_phreds,
        **kwargs,
    )
    consensus_seqs = poa_result["consensus_seqs"]
    if len(consensus_seqs) != 1:
        raise ValueError(
            f"expecting one consensus sequence, POA returned {len(consensus_seqs)}"
        )
    res = dict(consensus_seq=consensus_seqs[0])
    if return_phreds:
        res["consensus_phred"] = poa_result["consensus_phreds"][0]
    if names is not None:
        res["name"] = f"consensus_{_names_to_hash(names)}"
    return res


def get_consensus(df, use_phreds=True, return_phreds=True, return_name=True, **kwargs):
    seqs = df["read_seq"].to_list()
    rc = df["reverse_complement"].to_arrow()
    if use_phreds:
        phreds = df["read_phred"].to_arrow()
    else:
        phreds = None
    if return_name:
        names = df["name"].to_list()
    else:
        names = None
    return _get_consensus(
        seqs,
        rc,
        phreds=phreds,
        names=names,
        return_phreds=return_phreds,
        **kwargs,
    )


def get_consensus_group_by(
    df, use_phreds=True, return_phreds=True, return_name=True, **kwargs
):
    res = get_consensus(
        df,
        use_phreds=use_phreds,
        return_phreds=return_phreds,
        return_name=return_name,
        **kwargs,
    )
    # TODO: until polars accepts pyarrow/numpy arrays in pl.map_groups()
    return {k: list(v) if isinstance(v, np.ndarray) else v for k, v in res.items()}


# SEE: https://numba.discourse.group/t/feature-request-about-supporting-arrow-in-numba/1668/2
@numba.njit(nogil=True)
def align_phreds(seqs, phreds, gap_quality_method="mean"):
    num_seqs = len(seqs)
    if not num_seqs:
        return
    msa_length = len(seqs[0])
    aligned_phreds = np.empty((num_seqs, msa_length), dtype=np.int32)
    for seq_idx in range(num_seqs):
        last_nongap_phred = -1
        aligned_seq = seqs[seq_idx]
        aligned_phred = aligned_phreds[seq_idx]
        aligned_phred[:] = -1
        unaligned_phred = phreds[seq_idx]
        # necessary to avoid a NumbaTypeSafetyWarning arising from base_idx - offset
        offset = numba.types.int64(0)
        for base_idx in range(msa_length):
            base = aligned_seq[base_idx]
            if base == GAP_CHAR:
                offset += 1
                if last_nongap_phred != -1:
                    aligned_phred[base_idx] = last_nongap_phred
            else:
                phred = unaligned_phred[base_idx - offset]
                aligned_phred[base_idx] = phred
                last_nongap_phred = phred
        # let last_nongap_phred carry over
        # numba doesn't support reversed(range(msa_length))
        for base_idx in range(msa_length - 1, -1, -1):
            base = aligned_seq[base_idx]
            if base == GAP_CHAR:
                offset -= 1
                if last_nongap_phred != -1:
                    existing_aligned_phred = aligned_phred[base_idx]
                    if existing_aligned_phred == -1:
                        aligned_phred[base_idx] = last_nongap_phred
                    else:
                        if gap_quality_method == "min":
                            base_phred = min(last_nongap_phred, existing_aligned_phred)
                        elif gap_quality_method == "mean":
                            base_phred = (
                                last_nongap_phred + existing_aligned_phred
                            ) // 2
                        aligned_phred[base_idx] = base_phred
            else:
                base_phred = unaligned_phred[base_idx - offset]
                last_nongap_phred = base_phred
                aligned_phred[base_idx] = base_phred
    return aligned_phreds


def chars_to_str(ary, gaps=False):
    seq = ary.tobytes().decode()
    if not gaps:
        seq = seq.replace("-", "")
    return seq


def print_msa(seqs, phreds=None):
    if phreds is not None:
        aligned_phreds = align_phreds(seqs, phreds)
    msa_length = len(seqs[0])
    num_seqs = len(seqs)
    base_votes = []
    alphabet = set()
    max_votes = 0
    for base_idx in range(msa_length):
        votes = {}
        for seq_idx in range(num_seqs):
            base = seqs[seq_idx, base_idx]
            alphabet.add(base)
            if phreds is not None:
                vote = aligned_phreds[seq_idx, base_idx]
            else:
                vote = 1
            votes[base] = votes.get(base, 0) + vote
        max_votes = max(max_votes, max(votes.values()))
        base_votes.append(votes)
    digits = int(np.ceil(np.log10(max_votes)))
    alphabet = [*sorted(list(alphabet - set([GAP_CHAR]))), GAP_CHAR]
    print("".join(f"{' '*(digits)}{chr(b)}" for b in alphabet))
    for base_idx in range(msa_length):
        print(
            " "
            + " ".join(f"{{: >{digits}d}}" for _ in range(len(alphabet))).format(
                *tuple([base_votes[base_idx].get(b, 0) for b in alphabet])
            )
        )
