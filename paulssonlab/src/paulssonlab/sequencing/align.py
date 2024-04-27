from paulssonlab.sequencing.cigar import CigarOp, encode_cigar
from paulssonlab.util.core import pop_keys

try:
    import parasail
except:
    pass
try:
    from pywfa import WavefrontAligner
except:
    pass

PARASAIL_DEFAULTS = {"gap_opening": 10, "gap_extension": 1, "match": 1, "mismatch": -1}
PARASAIL_DEGENERATE_KWARGS = [
    "match",
    "mismatch",
    "degenerate_match",
    "degenerate_mismatch",
]
PARASAIL_ALGORITHMS = "nw,sg,sg_qb,sg_qe,sg_qx,sg_db,sg_de,sg_dx,sg_qb_de,sg_qe_db,sg_qb_db,sg_qe_de,sw".split(
    ","
)
PARASAIL_VECTORIZATION_STRATEGIES = ["striped", "scan", "diag"]
PARASAIL_SOLUTION_WIDTHS = ["8", "16", "32", "64", "sat"]

PYWFA_DEFAULTS = {}
PYWFA_CLASS_KWARGS = [
    "match",
    "gap_opening",
    "gap_extension",
    "gap_opening2",
    "gap_extension2",
    "scope",
    "span",
    "pattern_begin_free",
    "pattern_end_free",
    "text_begin_free",
    "text_end_free",
    "heuristic",
    "min_wavefront_length",
    "max_distance_threshold",
    "steps_between_cutoffs",
    "xdrop",
]
PYWFA_CALL_KWARGS = [
    "clip_cigar",
    "min_aligned_bases_left",
    "min_aligned_bases_right",
    "elide_mismatches",
    "supress_sequences",
]

DEGENERATE_BASES = {
    "R": "AG",
    "Y": "CT",
    "M": "AC",
    "K": "GT",
    "S": "CG",
    "W": "AT",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}


def degenerate_parasail_matrix(
    match,
    mismatch,
    degenerate_match=None,
    degenerate_mismatch=None,
    degenerate_bases=DEGENERATE_BASES,
):
    if degenerate_match is None:
        degenerate_match = match
    if degenerate_mismatch is None:
        degenerate_mismatch = mismatch
    bases = "ATCG" + "".join(degenerate_bases.keys())
    base_to_idx = {base: idx for idx, base in enumerate(bases)}
    if degenerate_match is None:
        degenerate_match = match
    if degenerate_mismatch is None:
        degenerate_mismatch = mismatch
    matrix = parasail.matrix_create(bases, match, mismatch)
    for deg_base, matching_bases in degenerate_bases.items():
        idx = base_to_idx[deg_base]
        degenerate_match_idxs = [base_to_idx[base] for base in matching_bases]
        degenerate_mismatch_idxs = [
            base_to_idx[base] for base in set("ATCG") - set(matching_bases)
        ]
        for idx2 in degenerate_match_idxs:
            matrix[idx, idx2] = matrix[idx2, idx] = degenerate_match
        for idx2 in degenerate_mismatch_idxs:
            matrix[idx, idx2] = matrix[idx2, idx] = degenerate_mismatch
    alphabet_aliases = "".join(
        f"{base}{deg_base}{deg_base}{base}"
        for deg_base, matching_bases in degenerate_bases.items()
        for base in matching_bases
    )
    return matrix, alphabet_aliases


def _decode_parasail_cigar(cigar_seq):
    # SEE: https://github.com/jeffdaily/parasail/blob/600fb26151ff19899ee39a214972dcf2b9b11ed7/src/cigar.c#L105
    return [(CigarOp(i & 0xF), i >> 4) for i in cigar_seq]


def _decode_pywfa_cigar(cigar_tuples):
    # use = instead of WFA2-lib's use of M to indicate base matches
    # (WFA2-lib correctly uses X to indicate mismatches)
    op_mapping = {0: 7}
    return [(CigarOp(op_mapping.get(op, op)), length) for op, length in cigar_tuples]


def pairwise_align(
    query,
    ref,
    method="parasail",
    parasail_algorithm="sw",
    parasail_vectorization_strategy="striped",
    parasail_solution_width="sat",
    parasail_case_sensitive=False,
    alphabet_aliases=None,
    degenerate=False,
    cigar_as_string=False,
    upper=True,
    **kwargs,
):
    if upper:
        ref = ref.upper()
        query = query.upper()
    if method == "parasail":
        kwargs = {**PARASAIL_DEFAULTS, **kwargs}
        gap_extension = kwargs.pop("gap_extension")
        gap_opening = kwargs.pop("gap_opening")
        degenerate_kwargs, kwargs = pop_keys(kwargs, PARASAIL_DEGENERATE_KWARGS)
        if kwargs:
            raise ValueError(f"unexpected kwargs: {kwargs}")
        if parasail_algorithm not in PARASAIL_ALGORITHMS:
            raise ValueError(f"invalid parasail algorithm: {parasail_algorithm}")
        if parasail_vectorization_strategy not in PARASAIL_VECTORIZATION_STRATEGIES:
            raise ValueError(
                f"invalid parasail vectorization strategy: {parasail_vectorization_strategy}"
            )
        parasail_solution_width = str(parasail_solution_width)  # in case int is given
        if parasail_solution_width not in PARASAIL_SOLUTION_WIDTHS:
            raise ValueError(
                f"invalid parasail solution width: {parasail_solution_width}"
            )
        parasail_function = f"{parasail_algorithm}_trace_{parasail_vectorization_strategy}_{parasail_solution_width}"
        parasail_align_func = getattr(parasail, parasail_function)
        if degenerate:
            substitution_matrix, alphabet_aliases = degenerate_parasail_matrix(
                **degenerate_kwargs
            )
        else:
            substitution_matrix = parasail.matrix_create(
                "ATCG", degenerate_kwargs["match"], degenerate_kwargs["mismatch"]
            )
        res = parasail_align_func(
            query, ref, gap_opening, gap_extension, substitution_matrix
        )
        score = res.score
        cigar = _decode_parasail_cigar(
            res.get_cigar(
                case_sensitive=parasail_case_sensitive,
                alphabet_aliases=alphabet_aliases,
            ).seq
        )
    elif method == "pywfa":
        kwargs = {**PYWFA_DEFAULTS, **kwargs}
        class_kwargs, kwargs = pop_keys(kwargs, PYWFA_CLASS_KWARGS)
        call_kwargs, kwargs = pop_keys(kwargs, PYWFA_CALL_KWARGS)
        if kwargs:
            raise ValueError(f"unexpected kwargs: {kwargs}")
        wfa = WavefrontAligner(**class_kwargs)
        wfa(query, ref, **call_kwargs)
        if wfa.status != 0:
            raise RuntimeError(f"alignment failed, pywfa returned status {wfa.status}")
        # WFA2-lib fixes matches to be worth 0 and mismatches/gaps/extensions have negative score (penalty)
        # SEE: https://github.com/smarco/WFA2-lib/issues/20
        score = wfa.score
        cigar = _decode_pywfa_cigar(wfa.cigartuples)
    else:
        raise ValueError(f"invalid method {method}")
    if cigar_as_string:
        cigar = encode_cigar(cigar)
    return score, cigar
