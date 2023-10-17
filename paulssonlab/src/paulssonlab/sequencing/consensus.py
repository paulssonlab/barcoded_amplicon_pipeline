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
SPOA_DEFAULTS = {}


def prepare_reads(seqs, rcs, phreds=None):
    if isinstance(seqs, pa.Array):
        seqs = seqs.to_pylist()  # prefer seqs to be given as a python list
    if isinstance(rcs, pa.Array):
        rcs = ak.from_arrow(rcs)
    seqs_oriented = [
        reverse_complement(seq) if rc else seq for seq, rc in zip(seqs, rcs)
    ]
    if phreds is not None:
        phreds_oriented = ak.from_arrow(
            pa.array(
                [
                    phred.values.to_numpy()[::-1] if rc else phred.values.to_numpy()
                    for phred, rc in zip(phreds, rcs)
                ]
            )
        )
        return seqs_oriented, phreds_oriented
    else:
        return seqs_oriented


def msa(seqs, method="abpoa", **kwargs):
    if method == "abpoa":
        aligner = pyabpoa.msa_aligner(**{**ABPOA_DEFAULTS, **kwargs})
        res = aligner.msa(seqs, out_cons=False, out_msa=True)
        msa_seqs = res.msa_seq
    elif method == "spoa":
        _, msa_seqs = spoa.poa(seqs, **{**SPOA_DEFAULTS, **kwargs})
    else:
        raise ValueError(f"method must be one of: abpoa, spoa")
    msa_seqs = np.array(
        [np.frombuffer(seq.encode(), dtype=np.uint8) for seq in msa_seqs]
    )
    return msa_seqs


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


# SEE: https://numba.discourse.group/t/feature-request-about-supporting-arrow-in-numba/1668/2
@numba.njit(nogil=True)
def phred_weighted_consensus(seqs, phreds, gap_quality_method="mean"):
    num_seqs = len(seqs)
    if not num_seqs:
        return None, None, None, None
    msa_length = len(seqs[0])
    aligned_phred = np.empty(msa_length, dtype=np.int32)
    votes = [
        numba.typed.Dict.empty(key_type=numba.types.uint8, value_type=numba.types.int32)
        for _ in range(msa_length)
    ]
    for seq_idx in range(num_seqs):
        last_nongap_phred = -1
        aligned_seq = seqs[seq_idx]
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
                        # not necessary, since we don't use aligned_phred
                        # aligned_phred[base_idx] = last_nongap_phred
                        base_phred = last_nongap_phred
                    else:
                        if gap_quality_method == "min":
                            base_phred = min(last_nongap_phred, existing_aligned_phred)
                        elif gap_quality_method == "mean":
                            base_phred = (
                                last_nongap_phred + existing_aligned_phred
                            ) // 2
                        # not necessary, since we don't use aligned_phred
                        # aligned_phred[base_idx] = base_phred
            else:
                base_phred = unaligned_phred[base_idx - offset]
                last_nongap_phred = base_phred
                # not necessary, since we don't use aligned_phred
                # aligned_phred[base_idx] = base_phred
            votes[base_idx][base] = votes[base_idx].get(base, 0) + base_phred
    consensus = np.empty(msa_length, dtype=np.uint8)
    nonconsensus = np.full(msa_length, SPACE_CHAR, dtype=np.uint8)
    consensus_phred = np.empty(msa_length, dtype=np.int32)
    nonconsensus_phred = np.zeros(msa_length, dtype=np.int32)
    for base_idx in range(msa_length):
        sorted_votes = sorted(
            [(v, k) for k, v in votes[base_idx].items()], reverse=True
        )
        consensus[base_idx] = sorted_votes[0][1]
        consensus_phred[base_idx] = sorted_votes[0][0]
        if len(sorted_votes) >= 2:
            nonconsensus[base_idx] = sorted_votes[1][1]
            # TODO: numba doesn't support sum
            # nonconsensus_phred[base_idx] = sum(v[0] for v in sorted_votes[1:])
            p = 0
            for v in sorted_votes[1:]:
                p += v[0]
            nonconsensus_phred[base_idx] = p
    return consensus, consensus_phred, nonconsensus, nonconsensus_phred


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
