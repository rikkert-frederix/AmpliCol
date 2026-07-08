from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Mapping, Sequence

from .color_plan import GenericColorPlan, LCColorSector

NC = 3


@dataclass(frozen=True)
class ColorContractionEntry:
    left_group_id: int
    right_group_id: int
    weight_re: float
    weight_im: float = 0.0
    symmetry_factor: float = 1.0

    def to_json_dict(self) -> dict[str, object]:
        return {
            "left_group_id": self.left_group_id,
            "right_group_id": self.right_group_id,
            "weight": [self.weight_re, self.weight_im],
            "symmetry_factor": self.symmetry_factor,
        }


@dataclass(frozen=True)
class ColorContractionPlan:
    color_accuracy: str
    supported: bool
    reason: str | None
    group_count: int
    entries: tuple[ColorContractionEntry, ...]
    includes_color_factor: bool = True

    def to_json_dict(self) -> dict[str, object]:
        return {
            "kind": "pyamplicol-color-contraction-plan",
            "color_accuracy": self.color_accuracy,
            "supported": self.supported,
            "reason": self.reason,
            "group_count": self.group_count,
            "includes_color_factor": self.includes_color_factor,
            "entry_count": len(self.entries),
            "storage": "upper-triangular sparse metric over coherent amplitude groups",
            "entries": [entry.to_json_dict() for entry in self.entries],
        }


@dataclass(frozen=True)
class ColorGroupDescriptor:
    group_id: int
    helicity_key: tuple[object, ...]
    sector_id: int
    word: tuple[int, ...]
    helicity_weight: float


def build_color_contraction_plan(
    color_plan: GenericColorPlan,
    groups: Sequence[ColorGroupDescriptor],
) -> ColorContractionPlan | None:
    """Build the final colour contraction over coherent amplitude groups.

    LC artifacts keep the historical Rusticol diagonal reduction path.  NLC and
    full-colour artifacts attach an explicit sparse colour matrix whose entries
    include AmpliCol's colour factors, so no leading-colour scalar factor should
    be applied again at runtime.
    """

    accuracy = color_plan.color_accuracy
    if accuracy == "lc":
        return None
    if accuracy not in {"nlc", "full"}:
        return ColorContractionPlan(
            color_accuracy=accuracy,
            supported=False,
            reason=f"unknown colour accuracy {accuracy!r}",
            group_count=len(groups),
            entries=(),
        )
    quark_pairs = color_plan.process.quark_lines.quark_pair_count
    if quark_pairs > 2:
        return ColorContractionPlan(
            color_accuracy=accuracy,
            supported=False,
            reason=(
                f"{accuracy} colour is currently matched to Fortran AmpliCol "
                "only for zero, one, or two quark pairs"
            ),
            group_count=len(groups),
            entries=(),
        )
    descriptors = tuple(groups)
    entries: list[ColorContractionEntry] = []
    sector_by_id = {sector.id: sector for sector in color_plan.sectors}
    factor_cache: dict[tuple[int, int], float] = {}
    descriptors_by_helicity: dict[tuple[object, ...], list[ColorGroupDescriptor]] = {}
    for descriptor in descriptors:
        descriptors_by_helicity.setdefault(descriptor.helicity_key, []).append(
            descriptor
        )
    for helicity_descriptors in descriptors_by_helicity.values():
        for left_offset, left in enumerate(helicity_descriptors):
            left_sector = sector_by_id.get(left.sector_id)
            if left_sector is None:
                return ColorContractionPlan(
                    color_accuracy=accuracy,
                    supported=False,
                    reason=f"missing colour sector {left.sector_id}",
                    group_count=len(groups),
                    entries=(),
                )
            for right in helicity_descriptors[left_offset:]:
                right_sector = sector_by_id.get(right.sector_id)
                if right_sector is None:
                    return ColorContractionPlan(
                        color_accuracy=accuracy,
                        supported=False,
                        reason=f"missing colour sector {right.sector_id}",
                        group_count=len(groups),
                        entries=(),
                    )
                sector_pair = (left.sector_id, right.sector_id)
                weight = factor_cache.get(sector_pair)
                if weight is None:
                    weight = amplicol_color_factor(
                        color_plan,
                        left_sector,
                        right_sector,
                        accuracy=accuracy,
                        full_col_acc=20,
                    )
                    factor_cache[sector_pair] = weight
                if abs(weight) <= 0.0:
                    continue
                symmetry = 1.0 if left.group_id == right.group_id else 2.0
                helicity_weight = _common_helicity_weight(left, right)
                entries.append(
                    ColorContractionEntry(
                        left_group_id=left.group_id,
                        right_group_id=right.group_id,
                        weight_re=helicity_weight * weight,
                        symmetry_factor=symmetry,
                    )
                )
    return ColorContractionPlan(
        color_accuracy=accuracy,
        supported=True,
        reason=None,
        group_count=len(groups),
        entries=tuple(entries),
    )


def amplicol_color_factors(
    color_plan: GenericColorPlan,
    left: LCColorSector,
    right: LCColorSector,
    *,
    full_col_acc: int = 20,
) -> tuple[float, float, float]:
    """Return AmpliCol-convention (LC, NLC, full) colour factors."""

    n_quark_pairs = color_plan.process.quark_lines.quark_pair_count
    n_ord = len(_coloured_word(left))
    if len(_coloured_word(right)) != n_ord:
        return (0.0, 0.0, 0.0)
    if n_quark_pairs == 0:
        return _pure_gluon_color_factors(left, right, n_ord, full_col_acc)
    if n_quark_pairs == 1:
        return _one_quark_line_color_factors(left, right, n_ord)
    if n_quark_pairs == 2:
        return _two_quark_line_color_factors(color_plan, left, right, n_ord)
    return (0.0, 0.0, 0.0)


def amplicol_color_factor(
    color_plan: GenericColorPlan,
    left: LCColorSector,
    right: LCColorSector,
    *,
    accuracy: str,
    full_col_acc: int = 20,
) -> float:
    """Return only the requested AmpliCol-convention colour factor."""

    n_quark_pairs = color_plan.process.quark_lines.quark_pair_count
    n_ord = len(_coloured_word(left))
    if len(_coloured_word(right)) != n_ord:
        return 0.0
    if n_quark_pairs == 0:
        return _pure_gluon_color_factor(
            left,
            right,
            n_ord,
            accuracy=accuracy,
            full_col_acc=full_col_acc,
        )
    values = amplicol_color_factors(
        color_plan,
        left,
        right,
        full_col_acc=full_col_acc,
    )
    if accuracy == "lc":
        return values[0]
    if accuracy == "nlc":
        return values[1]
    if accuracy == "full":
        return values[2]
    return 0.0


def _pure_gluon_color_factors(
    left: LCColorSector,
    right: LCColorSector,
    n_ord: int,
    full_col_acc: int,
) -> tuple[float, float, float]:
    return (
        _pure_gluon_color_factor(
            left,
            right,
            n_ord,
            accuracy="lc",
            full_col_acc=full_col_acc,
        ),
        _pure_gluon_color_factor(
            left,
            right,
            n_ord,
            accuracy="nlc",
            full_col_acc=full_col_acc,
        ),
        _pure_gluon_color_factor(
            left,
            right,
            n_ord,
            accuracy="full",
            full_col_acc=full_col_acc,
        ),
    )


def _pure_gluon_color_factor(
    left: LCColorSector,
    right: LCColorSector,
    n_ord: int,
    *,
    accuracy: str,
    full_col_acc: int,
) -> float:
    iper = _coloured_word(left)
    jper = _coloured_word(right)
    if accuracy == "lc":
        return float(NC**n_ord) if iper == jper else 0.0
    if accuracy == "nlc":
        if iper == jper:
            return float(NC**n_ord - n_ord * NC ** (n_ord - 2))
        return float(_check_nlc(tuple(jper), tuple(iper)) * NC ** (n_ord - 2))
    if accuracy != "full":
        return 0.0
    full_terms = _simplify_trace_terms(
        ((Fraction(1), (tuple(iper), tuple(reversed(jper)))),)
    )
    return _eval_nc_terms(
        full_terms,
        min_power=max(n_ord - 2 * full_col_acc, 0),
    )


def _one_quark_line_color_factors(
    left: LCColorSector,
    right: LCColorSector,
    n_ord: int,
) -> tuple[float, float, float]:
    iper = _coloured_word(left)
    jper = _coloured_word(right)
    lc = float(NC ** (n_ord - 1)) if iper == jper else 0.0
    full = _eval_trace(
        (tuple((*iper[1:-1], *reversed(jper[1:-1]))),),
    )
    if iper == jper:
        nlc = full
    else:
        acc = _check_nlc_1qqbar(tuple(jper[1:-1]), tuple(iper[1:-1]))
        nlc = full if acc != 0 else 0.0
    return (lc, nlc, full)


def _two_quark_line_color_factors(
    color_plan: GenericColorPlan,
    left: LCColorSector,
    right: LCColorSector,
    n_ord: int,
) -> tuple[float, float, float]:
    reference_start = _two_line_reference_start(color_plan)
    iper = _rotate_to_reference_start(_coloured_word(left), reference_start)
    jper = _rotate_to_reference_start(_coloured_word(right), reference_start)
    reference = (
        _rotate_to_reference_start(_coloured_word(color_plan.sectors[0]), reference_start)
        if color_plan.sectors
        else iper
    )
    gi, ui = _two_line_gi_ui(color_plan, iper, reference)
    gj, uj = _two_line_gi_ui(color_plan, jper, reference)
    same_flavour = _has_same_flavour_quark_lines(color_plan)
    lc = 0.0
    if iper == jper:
        if ui == 1 and uj == 1:
            lc = float(NC ** (n_ord - 2))
        elif ui == 2 and uj == 2 and not same_flavour:
            lc = float(NC ** (n_ord - 4) * 9.0)
        elif ui == 2 and uj == 2 and same_flavour:
            lc = float(NC ** (n_ord - 2))
    full = _two_quark_line_full_factor(
        iper,
        jper,
        n_ord=n_ord,
        gi=gi,
        gj=gj,
        ui=ui,
        uj=uj,
        same_flavour=same_flavour,
    )
    nlc = 0.0
    if abs(full) > 0.0:
        iper_glu, jper_glu = _two_line_ordered_gluon_strings(
            color_plan,
            iper,
            jper,
        )
        iper_ord, jper_ord = _convert_two_line_gluon_strings(
            n_ord,
            iper_glu,
            jper_glu,
        )
        acc = _check_nlc_2qqbar_sf(
            n_ord,
            iper_ord,
            jper_ord,
            gi,
            gj,
            ui,
            uj,
        )
        if acc != 0:
            nlc = full
    return (lc, nlc, full)


def _two_quark_line_full_factor(
    iper: tuple[int, ...],
    jper: tuple[int, ...],
    *,
    n_ord: int,
    gi: int,
    gj: int,
    ui: int,
    uj: int,
    same_flavour: bool,
) -> float:
    if ui == uj:
        traces = (
            tuple((*iper[1 : 1 + gi], *reversed(jper[1 : 1 + gj]))),
            tuple((*iper[gi + 3 : n_ord - 1], *reversed(jper[gj + 3 : n_ord - 1]))),
        )
        coeff = Fraction(1)
    elif (ui, uj) in {(1, 2), (2, 1)}:
        traces = (
            tuple(
                (
                    *iper[1 : 1 + gi],
                    *reversed(jper[gj + 3 : n_ord - 1]),
                    *iper[gi + 3 : n_ord - 1],
                    *reversed(jper[1 : 1 + gj]),
                )
            ),
        )
        coeff = Fraction(-1)
    else:
        return 0.0
    return _eval_trace(traces, coeff=coeff)


def _coloured_word(sector: LCColorSector) -> tuple[int, ...]:
    if sector.kind == "single-trace":
        return tuple(sector.trace_labels)
    if sector.kind == "open-lines":
        return tuple(sector.word_labels or sector.color_words[0])
    return ()


def _two_line_reference_start(color_plan: GenericColorPlan) -> int | None:
    if color_plan.process.quark_labels:
        return int(color_plan.process.quark_labels[0])
    return None


def _rotate_to_reference_start(
    word: tuple[int, ...],
    reference_start: int | None,
) -> tuple[int, ...]:
    if reference_start is None or reference_start not in word:
        return word
    offset = word.index(reference_start)
    return word[offset:] + word[:offset]


def _two_line_gi_ui(
    color_plan: GenericColorPlan,
    word: tuple[int, ...],
    reference_word: tuple[int, ...],
) -> tuple[int, int]:
    quarkish = set(color_plan.process.quark_labels) | set(
        color_plan.process.antiquark_labels
    )
    gi = 0
    for position in range(1, max(len(word) - 1, 1)):
        if word[position] in quarkish:
            gi = position - 1
            break
    ui = 1
    if len(word) >= 2 and len(reference_word) >= 2:
        same_ends = word[0] == reference_word[0] and word[-1] == reference_word[-1]
        opposite_ends = word[0] != reference_word[0] and word[-1] != reference_word[-1]
        ui = 1 if same_ends or opposite_ends else 2
    return gi, ui


def _two_line_ordered_gluon_strings(
    color_plan: GenericColorPlan,
    iper: tuple[int, ...],
    jper: tuple[int, ...],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    quarkish = set(color_plan.process.quark_labels) | set(
        color_plan.process.antiquark_labels
    )
    return (
        tuple(label for label in iper if label not in quarkish),
        tuple(label for label in jper if label not in quarkish),
    )


def _convert_two_line_gluon_strings(
    n_ord: int,
    iper: tuple[int, ...],
    jper: tuple[int, ...],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Mirror AmpliCol's convert_gluon_string helper.

    The NLC topology tests depend only on the relative order of gluons.  The
    Fortran routine maps the largest external label to 1, the next largest to
    2, and so on, separately for the row and column strings.
    """

    expected = max(n_ord - 4, 0)
    if len(iper) != expected or len(jper) != expected:
        return iper, jper

    def convert(word: tuple[int, ...]) -> tuple[int, ...]:
        ordered = sorted(word, reverse=True)
        rank = {label: index + 1 for index, label in enumerate(ordered)}
        return tuple(rank[label] for label in word)

    return convert(iper), convert(jper)


def _check_nlc_2qqbar_sf(
    n_ord: int,
    iper: tuple[int, ...],
    jper: tuple[int, ...],
    ri: int,
    rj: int,
    ii: int,
    jj: int,
) -> int:
    """Port of AmpliCol's check_NLC_2qqbar_SF topology filter.

    It returns 99 for LC-like entries, +/-1 for first subleading entries, and
    0 for NNLC-or-beyond entries.  The final NLC coefficient is still the full
    trace overlap; this function only decides whether that overlap belongs in
    the NLC-expanded matrix.
    """

    m = max(n_ord - 4, 0)
    if len(iper) != m or len(jper) != m:
        return 0
    if (ii, jj) in {(1, 2), (2, 1)}:
        temp = (
            iper[:ri]
            + tuple(reversed(jper[rj:]))
            + iper[ri:]
            + tuple(reversed(jper[:rj]))
        )
        return -1 if _nlc_pairing_is_planar(temp, m) else 0
    if (ii, jj) not in {(1, 1), (2, 2)}:
        return 0
    if ri == rj and iper == jper:
        return 99

    aa = iper[:ri]
    bb = iper[ri:]
    cc = jper[:rj]
    dd = jper[rj:]
    temp2 = aa + tuple(reversed(cc))
    temp3 = bb + tuple(reversed(dd))
    disjoint = not (set(temp2) & set(temp3))
    if disjoint:
        if aa and (bb == dd or not bb):
            return _check_nlc_subword(cc, aa, threshold=m - 4)
        if bb and (aa == cc or not aa):
            return _check_nlc_subword(dd, bb, threshold=m - rj - 4)
        return 0
    return _check_nlc_2qqbar_overlap(aa, bb, cc, dd, m, ri, rj)


def _nlc_pairing_is_planar(word: tuple[int, ...], m: int) -> bool:
    if m == 0:
        return True
    positions: dict[int, list[int]] = {label: [] for label in range(1, m + 1)}
    for offset, label in enumerate(word):
        if label in positions:
            positions[label].append(offset)
    if any(len(item) != 2 for item in positions.values()):
        return False
    for first, pair in positions.items():
        if abs(pair[0] - pair[1]) % 2 != 1:
            return False
        for second in range(first + 1, m + 1):
            if not _intervals_disjoint_or_nested(tuple(pair), tuple(positions[second])):
                return False
    return True


def _intervals_disjoint_or_nested(
    first: tuple[int, int],
    second: tuple[int, int],
) -> bool:
    a1, a2 = first
    b1, b2 = second
    return (
        (a1 < b1 and a2 < b1)
        or (a1 > b2 and a2 > b2)
        or (a1 > b1 and a2 < b2)
        or (b1 > a1 and b2 < a2)
    )


def _check_nlc_subword(
    candidate: tuple[int, ...],
    target: tuple[int, ...],
    *,
    threshold: int,
) -> int:
    if candidate == target:
        return 99
    if not candidate or len(candidate) != len(target):
        return 0
    sign = _check_nlc(candidate, target)
    if sign == 0:
        return 0
    if sign < 0:
        return sign
    return sign


def _check_nlc_2qqbar_overlap(
    aa: tuple[int, ...],
    bb: tuple[int, ...],
    cc: tuple[int, ...],
    dd: tuple[int, ...],
    m: int,
    ri: int,
    rj: int,
) -> int:
    """Remaining overlapping-sets branch of check_NLC_2qqbar_SF."""

    temp2 = aa + tuple(reversed(cc))
    temp3 = bb + tuple(reversed(dd))
    common = set(temp2) & set(temp3)
    if not common:
        return 0
    skipped = None
    ind_i = ind_j = -1
    for i, label in enumerate(temp2):
        if label not in common:
            continue
        ind_i = i
        ind_j = temp3.index(label)
        skipped = label
        break
    if skipped is None:
        return 0
    if len(temp2) == 1 or len(temp3) == 1:
        return 0

    perm: tuple[int, ...]
    if ind_i < ri and ind_j >= len(bb):
        # Common generator in the A-D pair.
        perm = (
            temp2[ind_i + 1 : ri]
            + temp2[ri : ri + rj]
            + temp2[:ind_i]
            + temp3[ind_j + 1 :]
            + temp3[: len(bb)]
            + temp3[len(bb) : ind_j]
        )
        itemp4 = temp2[:ind_i] + temp2[ind_i + 1 : ri]
        itemp5 = tuple(reversed(temp3[ind_j + 1 :])) + tuple(
            reversed(temp3[len(bb) : ind_j])
        )
        if rj == ri - 1 and itemp4 == cc:
            return 0
        if rj + 1 == ri and itemp5 == bb:
            return 0
    elif ind_i >= ri and ind_j < len(bb):
        # Common generator in the B-C pair.
        perm = (
            temp2[ind_i + 1 :]
            + temp2[:ri]
            + temp2[ri:ind_i]
            + temp3[ind_j + 1 : len(bb)]
            + temp3[len(bb) :]
            + temp3[:ind_j]
        )
        itemp6 = tuple(reversed(temp2[ind_i + 1 :])) + tuple(
            reversed(temp2[ri:ind_i])
        )
        itemp7 = temp3[:ind_j] + temp3[ind_j + 1 : len(bb)]
        if rj - 1 == ri and itemp6 == aa:
            return 0
        if rj == ri + 1 and itemp7 == dd:
            return 0
    else:
        return 0

    positions: dict[int, list[int]] = {label: [] for label in range(1, m + 1)}
    for offset, label in enumerate(perm):
        if label == skipped:
            continue
        if label in positions:
            positions[label].append(offset)
    for label, pair in positions.items():
        if label == skipped:
            continue
        if len(pair) != 2 or abs(pair[0] - pair[1]) % 2 != 1:
            return 0
    for first in range(1, m + 1):
        if first == skipped:
            continue
        for second in range(first + 1, m + 1):
            if second == skipped:
                continue
            if not _intervals_disjoint_or_nested(
                tuple(positions[first]),
                tuple(positions[second]),
            ):
                return 0
    return 1


def _has_same_flavour_quark_lines(color_plan: GenericColorPlan) -> bool:
    pdg_by_label: dict[int, int] = {
        leg.label: abs(int(leg.outgoing_pdg))
        for leg in color_plan.process.legs
        if leg.outgoing_pdg is not None
    }
    quark_flavours = [
        pdg_by_label[label]
        for label in color_plan.process.quark_labels
        if label in pdg_by_label
    ]
    return len(set(quark_flavours)) < len(quark_flavours)


def _common_helicity_weight(
    left: ColorGroupDescriptor,
    right: ColorGroupDescriptor,
) -> float:
    if abs(left.helicity_weight - right.helicity_weight) > 1.0e-12:
        return 0.5 * (left.helicity_weight + right.helicity_weight)
    return left.helicity_weight


def _eval_trace(
    traces: tuple[tuple[int, ...], ...],
    *,
    coeff: Fraction = Fraction(1),
) -> float:
    terms = _simplify_trace_terms(((coeff, traces),))
    return _eval_nc_terms(terms)


def _simplify_trace_terms(
    initial: Iterable[tuple[Fraction, tuple[tuple[int, ...], ...]]],
) -> Mapping[int, Fraction]:
    terms: list[tuple[Fraction, tuple[tuple[int, ...], ...]]] = [
        (coeff, tuple(tuple(trace) for trace in traces))
        for coeff, traces in initial
        if coeff
    ]
    guard = 0
    while True:
        guard += 1
        if guard > 100000:
            raise RuntimeError("colour trace simplification did not converge")
        terms = _simplify_tr1_tr0(terms)
        if not terms:
            return {}
        if all(not traces for _coeff, traces in terms):
            result: dict[int, Fraction] = {}
            for coeff, _traces in terms:
                result[0] = result.get(0, Fraction(0)) + coeff
            return result
        next_terms, changed = _trace_pair_simplify(terms)
        terms = _simplify_tr1_tr0(next_terms)
        next_terms, changed_dup = _trace_duplicate_simplify(terms)
        terms = next_terms
        if not changed and not changed_dup and any(traces for _coeff, traces in terms):
            raise RuntimeError(f"cannot simplify colour traces: {terms[:3]}")


def _simplify_tr1_tr0(
    terms: list[tuple[Fraction, tuple[tuple[int, ...], ...]]],
) -> list[tuple[Fraction, tuple[tuple[int, ...], ...]]]:
    result: list[tuple[Fraction, tuple[tuple[int, ...], ...]]] = []
    for coeff, traces in terms:
        if any(len(trace) == 1 for trace in traces):
            continue
        power = sum(1 for trace in traces if len(trace) == 0)
        remaining = tuple(trace for trace in traces if len(trace) != 0)
        result.append((coeff * (NC**power), remaining))
    return [(coeff, traces) for coeff, traces in result if coeff]


def _trace_pair_simplify(
    terms: list[tuple[Fraction, tuple[tuple[int, ...], ...]]],
) -> tuple[list[tuple[Fraction, tuple[tuple[int, ...], ...]]], bool]:
    result: list[tuple[Fraction, tuple[tuple[int, ...], ...]]] = []
    changed = False
    for coeff, traces in terms:
        replaced = False
        for first_index, first in enumerate(traces):
            for second_index in range(first_index + 1, len(traces)):
                second = traces[second_index]
                found = _first_common_position(first, second)
                if found is None:
                    continue
                first_pos, second_pos = found
                combined = (
                    first[:first_pos]
                    + second[second_pos + 1 :]
                    + second[:second_pos]
                    + first[first_pos + 1 :]
                )
                rest = tuple(
                    trace
                    for index, trace in enumerate(traces)
                    if index not in {first_index, second_index}
                )
                result.append((coeff, (combined, *rest)))
                reduced_first = first[:first_pos] + first[first_pos + 1 :]
                reduced_second = second[:second_pos] + second[second_pos + 1 :]
                result.append(
                    (
                        -coeff / NC,
                        tuple(
                            trace
                            if index not in {first_index, second_index}
                            else (reduced_first if index == first_index else reduced_second)
                            for index, trace in enumerate(traces)
                        ),
                    )
                )
                changed = True
                replaced = True
                break
            if replaced:
                break
        if not replaced:
            result.append((coeff, traces))
    return result, changed


def _trace_duplicate_simplify(
    terms: list[tuple[Fraction, tuple[tuple[int, ...], ...]]],
) -> tuple[list[tuple[Fraction, tuple[tuple[int, ...], ...]]], bool]:
    result: list[tuple[Fraction, tuple[tuple[int, ...], ...]]] = []
    changed = False
    for coeff, traces in terms:
        replaced = False
        for trace_index, trace in enumerate(traces):
            positions = _first_duplicate_positions(trace)
            if positions is None:
                continue
            first_pos, second_pos = positions
            a = trace[:first_pos]
            b = trace[first_pos + 1 : second_pos]
            c = trace[second_pos + 1 :]
            rest = tuple(
                item for index, item in enumerate(traces) if index != trace_index
            )
            result.append((coeff, (a + c, b, *rest)))
            result.append((-coeff / NC, (a + b + c, *rest)))
            changed = True
            replaced = True
            break
        if not replaced:
            result.append((coeff, traces))
    return result, changed


def _first_common_position(
    left: Sequence[int],
    right: Sequence[int],
) -> tuple[int, int] | None:
    for left_index, label in enumerate(left):
        for right_index, other in enumerate(right):
            if label == other:
                return left_index, right_index
    return None


def _first_duplicate_positions(trace: Sequence[int]) -> tuple[int, int] | None:
    seen: dict[int, int] = {}
    for index, label in enumerate(trace):
        if label in seen:
            return seen[label], index
        seen[label] = index
    return None


def _eval_nc_terms(
    terms: Mapping[int, Fraction],
    *,
    min_power: int | None = None,
) -> float:
    value = 0.0
    for power, coeff in terms.items():
        if min_power is not None and power < min_power:
            continue
        value += float(coeff) * (float(NC) ** power)
    return value


def _check_nlc(jper: tuple[int, ...], iper: tuple[int, ...]) -> int:
    n = len(iper)
    if n == 0:
        return 99
    i1 = next((i for i in range(n) if jper[i] != iper[i]), n)
    if i1 == n:
        return 99
    try:
        i2 = next(i for i in range(i1 + 1, n) if jper[i] == iper[i1])
    except StopIteration:
        return 0
    offset = 0
    while i2 + offset < n and i1 + offset < n and jper[i2 + offset] == iper[i1 + offset]:
        offset += 1
    i3 = i2 + offset - 1
    i4 = i1 + offset - 1
    if i4 + 1 >= n:
        return 0
    try:
        i5 = next(i for i in range(i1, n) if jper[i] == iper[i4 + 1])
    except StopIteration:
        return 0
    if i5 > i3:
        return 0
    sign = 1
    if i1 > i5 - 1:
        left_len = i4 - i1
        right_len = i2 - 1 - i5
        if left_len == 0 and right_len == 0:
            sign = -1
        elif n > 0 and (
            (left_len == 0 and right_len == n - 3)
            or (left_len == n - 3 and right_len == 0)
        ):
            sign = -1
        elif left_len == 0 or right_len == 0:
            return 0
        elif left_len + right_len <= max(n - 4, 0):
            pass
        else:
            return 0
    itemp = (
        jper[:i1]
        + jper[i2 : i3 + 1]
        + jper[i5:i2]
        + jper[i1:i5]
        + jper[i3 + 1 :]
    )
    return sign if itemp == iper else 0


def _check_nlc_1qqbar(jper: tuple[int, ...], iper: tuple[int, ...]) -> int:
    """Port AmpliCol's check_NLC_1qqbar for one open quark line."""

    n = len(iper)
    if len(jper) != n:
        return 0
    i1 = next((i for i in range(n) if jper[i] != iper[i]), n)
    if i1 >= n:
        return 99
    try:
        i2 = next(i for i in range(i1 + 1, n) if jper[i] == iper[i1])
    except StopIteration:
        return 0
    offset = 0
    while i2 + offset < n and i1 + offset < n and jper[i2 + offset] == iper[i1 + offset]:
        offset += 1
    i3 = i2 + offset - 1
    i4 = i1 + offset - 1
    if i4 + 1 >= n:
        return 0
    try:
        i5 = next(i for i in range(i1, n) if jper[i] == iper[i4 + 1])
    except StopIteration:
        return 0
    if i5 > i3:
        return 0
    sign = 1
    if i1 > i5 - 1:
        left_len = i4 - i1
        right_len = i2 - 1 - i5
        if left_len == 0 and right_len == 0:
            sign = -1
        elif left_len == 0 or right_len == 0:
            return 0
        elif left_len + right_len <= n - 4:
            pass
        else:
            sign = 1
    itemp = (
        jper[:i1]
        + jper[i2 : i3 + 1]
        + jper[i5:i2]
        + jper[i1:i5]
        + jper[i3 + 1 :]
    )
    return sign if itemp == iper else 0
