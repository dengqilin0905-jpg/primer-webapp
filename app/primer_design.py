import re
from typing import Any

from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

SNAPGENE_PRESET = {
    "Na": 50.0,
    "Mg": 1.5,
    "dNTPs": 0.2,
    "primer_conc": 250.0,
}
STOP_CODONS = {"TAA", "TAG", "TGA"}


def normalize_seq(s: str) -> str:
    s = (s or "").upper()
    s = re.sub(r"[^ACGT]", "", s)
    return s


def revcomp(s: str) -> str:
    return str(Seq(s).reverse_complement())


def trim_stop_codon(cds: str) -> tuple[str, bool]:
    cds = normalize_seq(cds)
    if len(cds) < 3:
        return cds, False
    stop = cds[-3:]
    if stop in STOP_CODONS:
        return cds[:-3], True
    return cds, False


def tm_nn(seq: str, Na=50.0, Mg=1.5, dNTPs=0.2, primer_conc=250.0):
    return mt.Tm_NN(seq, Na=Na, Mg=Mg, dNTPs=dNTPs, dnac1=primer_conc, dnac2=primer_conc)


def _contains_forbidden_homopolymer(s: str) -> bool:
    ss = (s or "").upper()
    return ("AAAAA" in ss) or ("TTTTT" in ss)


def _in_range(x: float, lo: float, hi: float) -> bool:
    return lo <= x <= hi


def _distance_to_range(x: float, lo: float, hi: float) -> float:
    if x < lo:
        return lo - x
    if x > hi:
        return x - hi
    return 0.0


def _simple_structure_penalty(primer: str) -> float:
    # Lightweight heuristic only: keep as soft ranking, not a hard rejection.
    p = (primer or "").upper()
    if not p:
        return 10.0

    penalty = 0.0
    if p[-1] in {"A", "T"}:
        penalty += 0.6
    if len(p) >= 2 and p[-2:] in {"AA", "TT"}:
        penalty += 0.3

    if len(p) >= 5:
        tail5 = p[-5:]
        rc_tail5 = revcomp(tail5)
        if rc_tail5 in p[:-5]:
            penalty += 1.0
    return penalty


def _len_pref_distance(length: int, pref_lo: int, pref_hi: int) -> float:
    if length < pref_lo:
        return float(pref_lo - length)
    if length > pref_hi:
        return float(length - pref_hi)
    return 0.0


def _candidate_score(
    primer: str,
    tm: float,
    gc: float,
    target_tm: float,
    gc_min: float,
    gc_max: float,
    pref_len_lo: int,
    pref_len_hi: int,
    tm_priority_weight: float,
) -> tuple[float, float, float, float, float]:
    tm_dev = abs(tm - target_tm)
    gc_dev = _distance_to_range(gc, gc_min, gc_max)
    struct_pen = _simple_structure_penalty(primer)
    len_dev = _len_pref_distance(len(primer), pref_len_lo, pref_len_hi)

    # Priority: Tm > GC > structure > length (length is soft preference only).
    score = (
        tm_priority_weight * tm_dev
        + 0.9 * gc_dev
        + 0.4 * struct_pen
        + 0.08 * len_dev
    )
    return score, tm_dev, gc_dev, struct_pen, len_dev


def _empty_stats() -> dict[str, int]:
    return {"total": 0, "pass": 0, "tm_fail": 0, "gc_fail": 0, "hp_fail": 0}


def _merge_stats(dst: dict[str, int], src: dict[str, int]):
    for k in dst.keys():
        dst[k] += int(src.get(k, 0))


def _resolve_tm_profile(
    tm_profile: str,
    Na: float,
    Mg: float,
    dNTPs: float,
    primer_conc: float,
) -> tuple[str, dict[str, float]]:
    norm = (tm_profile or "custom").strip().lower()
    if norm in {"snapgene", "snapgene_preset"}:
        return "snapgene_preset", dict(SNAPGENE_PRESET)
    return "custom", {
        "Na": float(Na),
        "Mg": float(Mg),
        "dNTPs": float(dNTPs),
        "primer_conc": float(primer_conc),
    }


def _iter_forward_candidates_at_start(
    seq: str,
    start: int,
    min_len: int,
    max_len: int,
    tm_min: float,
    tm_max: float,
    target_tm: float,
    gc_min: float,
    gc_max: float,
    pref_len_lo: int,
    pref_len_hi: int,
    tm_priority_weight: float,
    forbid_homopolymer: bool,
    tm_params: dict[str, float],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    n = len(seq)
    stats = _empty_stats()
    out: list[dict[str, Any]] = []
    for L in range(min_len, max_len + 1):
        end = start + L
        if end > n:
            break
        bind = seq[start:end]
        stats["total"] += 1
        t = tm_nn(bind, **tm_params)
        if not _in_range(t, tm_min, tm_max):
            stats["tm_fail"] += 1
            continue
        gc = float(gc_fraction(bind) * 100.0)
        if not _in_range(gc, gc_min, gc_max):
            stats["gc_fail"] += 1
            continue
        if forbid_homopolymer and _contains_forbidden_homopolymer(bind):
            stats["hp_fail"] += 1
            continue
        score, tm_dev, gc_dev, struct_pen, len_dev = _candidate_score(
            primer=bind,
            tm=t,
            gc=gc,
            target_tm=target_tm,
            gc_min=gc_min,
            gc_max=gc_max,
            pref_len_lo=pref_len_lo,
            pref_len_hi=pref_len_hi,
            tm_priority_weight=tm_priority_weight,
        )
        stats["pass"] += 1
        out.append(
            {
                "bind_start": start,
                "bind_end": end,
                "length": L,
                "primer_seq": bind,
                "tm": float(t),
                "gc": gc,
                "score": float(score),
                "tm_dev": float(tm_dev),
                "gc_dev": float(gc_dev),
                "struct_pen": float(struct_pen),
                "len_dev": float(len_dev),
            }
        )
    return out, stats


def _iter_reverse_candidates_at_end(
    seq: str,
    end: int,
    min_len: int,
    max_len: int,
    tm_min: float,
    tm_max: float,
    target_tm: float,
    gc_min: float,
    gc_max: float,
    pref_len_lo: int,
    pref_len_hi: int,
    tm_priority_weight: float,
    forbid_homopolymer: bool,
    tm_params: dict[str, float],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    stats = _empty_stats()
    out: list[dict[str, Any]] = []
    for L in range(min_len, max_len + 1):
        start = end - L
        if start < 0:
            continue
        bind = seq[start:end]
        rc = revcomp(bind)
        stats["total"] += 1
        t = tm_nn(rc, **tm_params)
        if not _in_range(t, tm_min, tm_max):
            stats["tm_fail"] += 1
            continue
        gc = float(gc_fraction(rc) * 100.0)
        if not _in_range(gc, gc_min, gc_max):
            stats["gc_fail"] += 1
            continue
        if forbid_homopolymer and _contains_forbidden_homopolymer(rc):
            stats["hp_fail"] += 1
            continue
        score, tm_dev, gc_dev, struct_pen, len_dev = _candidate_score(
            primer=rc,
            tm=t,
            gc=gc,
            target_tm=target_tm,
            gc_min=gc_min,
            gc_max=gc_max,
            pref_len_lo=pref_len_lo,
            pref_len_hi=pref_len_hi,
            tm_priority_weight=tm_priority_weight,
        )
        stats["pass"] += 1
        out.append(
            {
                "bind_start": start,
                "bind_end": end,
                "length": L,
                "primer_seq": rc,
                "tm": float(t),
                "gc": gc,
                "score": float(score),
                "tm_dev": float(tm_dev),
                "gc_dev": float(gc_dev),
                "struct_pen": float(struct_pen),
                "len_dev": float(len_dev),
            }
        )
    return out, stats


def _build_failure_note(prefix: str, stats: dict[str, int], tm_min: float, tm_max: float, gc_min: float, gc_max: float) -> str:
    return (
        f"{prefix} no candidate passing Tm[{tm_min:.1f},{tm_max:.1f}] and GC[{gc_min:.1f},{gc_max:.1f}]"
        f" (checked={stats['total']}, tm_fail={stats['tm_fail']}, gc_fail={stats['gc_fail']}, hp_fail={stats['hp_fail']})"
    )


def _pick_best_by_score(candidates: list[dict[str, Any]]) -> dict[str, Any] | None:
    if not candidates:
        return None
    return min(
        candidates,
        key=lambda c: (
            c["score"],
            c["tm_dev"],
            c["gc_dev"],
            c["struct_pen"],
            c["len_dev"],
            -c["length"],
        ),
    )


def _pair_key(
    fw: dict[str, Any],
    rv: dict[str, Any],
    tm_pair_delta_max: float,
    prefer_longest: bool,
) -> tuple[Any, ...]:
    delta = abs(float(fw["tm"]) - float(rv["tm"]))
    delta_pen = max(0.0, delta - tm_pair_delta_max)
    quality = float(fw["score"]) + float(rv["score"]) + 0.25 * delta + 2.5 * delta_pen
    if prefer_longest:
        amp = int(rv["bind_end"]) - int(fw["bind_start"])
        return (-amp, quality, delta, int(fw["bind_start"]), -int(rv["bind_end"]))
    return (quality, delta, -(int(fw["length"]) + int(rv["length"])))


def _select_pair(
    forwards: list[dict[str, Any]],
    reverses: list[dict[str, Any]],
    tm_pair_delta_max: float,
    tm_delta_mode: str,
    prefer_longest: bool,
) -> dict[str, Any] | None:
    if not forwards or not reverses:
        return None

    mode = (tm_delta_mode or "soft").strip().lower()
    if mode not in {"soft", "hard"}:
        mode = "soft"

    best_ok: dict[str, Any] | None = None
    best_any: dict[str, Any] | None = None

    for fw in forwards:
        for rv in reverses:
            if int(fw["bind_end"]) > int(rv["bind_start"]):
                continue
            delta = abs(float(fw["tm"]) - float(rv["tm"]))
            ok = delta <= tm_pair_delta_max + 1e-9
            if mode == "hard" and not ok:
                continue

            key = _pair_key(fw, rv, tm_pair_delta_max, prefer_longest)
            item = {
                "fw": fw,
                "rv": rv,
                "tm_delta": float(delta),
                "amplicon_len": int(rv["bind_end"]) - int(fw["bind_start"]),
                "key": key,
                "delta_relaxed": bool(not ok),
            }
            if best_any is None or item["key"] < best_any["key"]:
                best_any = item
            if ok and (best_ok is None or item["key"] < best_ok["key"]):
                best_ok = item

    if mode == "hard":
        return best_ok
    if best_ok is not None:
        best_ok["delta_relaxed"] = False
        return best_ok
    return best_any


def _build_cds_attempts(tm_min: float, tm_max: float, max_len: int) -> list[dict[str, Any]]:
    attempts: list[dict[str, Any]] = [{"tier": "tier0", "tm_min": float(tm_min), "tm_max": float(tm_max), "max_len": int(max_len)}]

    current = int(max_len)
    while current < 36:
        current = min(36, current + 2)
        attempts.append(
            {
                "tier": f"tier1_len{current}",
                "tm_min": float(tm_min),
                "tm_max": float(tm_max),
                "max_len": int(current),
            }
        )

    tm2_min = min(float(tm_min), 57.0)
    tm2_max = max(float(tm_max), 63.0)
    attempts.append({"tier": "tier2_tm57_63", "tm_min": tm2_min, "tm_max": tm2_max, "max_len": max(36, int(max_len))})

    tm3_min = min(tm2_min, 56.0)
    tm3_max = max(tm2_max, 64.0)
    attempts.append({"tier": "tier3_tm56_64", "tm_min": tm3_min, "tm_max": tm3_max, "max_len": max(36, int(max_len))})

    # Deduplicate attempts when user already provides broad ranges.
    uniq = []
    seen = set()
    for a in attempts:
        key = (round(a["tm_min"], 4), round(a["tm_max"], 4), int(a["max_len"]))
        if key in seen:
            continue
        seen.add(key)
        uniq.append(a)
    return uniq


def _frame_check_messages(
    seq: str,
    seq_type_norm: str,
    construct_type: str,
    remove_stop: bool,
    restriction_site_f: str,
    restriction_site_r: str,
    junction_extra_bases_f: int,
    junction_extra_bases_r: int,
) -> list[str]:
    msgs: list[str] = []
    ct = (construct_type or "CDS_no_tag").strip()
    ct_lower = ct.lower()
    if seq_type_norm != "CDS" or ct_lower == "promoter_or_non_coding":
        return msgs

    if len(seq) % 3 == 0:
        msgs.append("OK: CDS length is in-frame.")
    else:
        msgs.append("Warning: insert length is not a multiple of 3.")

    codons = [seq[i : i + 3] for i in range(0, max(0, len(seq) - 2), 3)]
    terminal_is_stop = bool(codons) and codons[-1] in STOP_CODONS
    internal = [c for c in codons[:-1] if c in STOP_CODONS] if terminal_is_stop else [c for c in codons if c in STOP_CODONS]
    if internal:
        msgs.append("Warning: internal stop codon detected.")

    if terminal_is_stop and not remove_stop:
        msgs.append("Warning: terminal stop retained; C-terminal fusion may fail.")
        if ct_lower == "c_terminal_tag":
            msgs.append("Warning: terminal stop retained; C-terminal tag will not be translated.")

    site_f = normalize_seq(restriction_site_f or "")
    site_r = normalize_seq(restriction_site_r or "")
    j5 = max(0, int(junction_extra_bases_f or 0))
    j3 = max(0, int(junction_extra_bases_r or 0))
    add5 = len(site_f) + j5
    add3 = len(site_r) + j3

    if add5 % 3 != 0:
        msgs.append(f"Warning: 5' junction adds {add5} bp, causing +{add5 % 3} frameshift.")
    if (len(seq) + add5 + add3) % 3 != 0:
        msgs.append("Warning: junction + insert is not in-frame.")
    elif add5 == 0 and add3 == 0 and len(seq) % 3 == 0:
        msgs.append("OK: junction + insert remains in-frame.")
    return msgs


def design_primers_for_gene(
    gene_name: str,
    seq: str,
    homo_f: str,
    homo_r: str,
    remove_stop: bool = True,
    min_len: int = 18,
    max_len: int = 30,
    tm_min: float = 55.0,
    tm_max: float = 65.0,
    target_tm: float = 60.0,
    gc_min: float = 40.0,
    gc_max: float = 60.0,
    Na: float = 50.0,
    Mg: float = 1.5,
    dNTPs: float = 0.2,
    primer_conc: float = 250.0,
    seq_type: str = "CDS",
    length_mode: str = "adaptive_soft",
    cds_len_pref: tuple[int, int] = (18, 21),
    non_cds_len_pref: tuple[int, int] = (20, 26),
    tm_priority_weight: float = 2.5,
    non_cds_forbid_homopolymer: bool = True,
    tm_pair_delta_max: float = 2.0,
    tm_delta_mode: str = "soft",
    tm_profile: str = "custom",
    cds_recovery_enable: bool = True,
    enable_frame_check: bool = False,
    construct_type: str = "CDS_no_tag",
    restriction_site_f: str = "",
    restriction_site_r: str = "",
    junction_extra_bases_f: int = 0,
    junction_extra_bases_r: int = 0,
):
    seq = normalize_seq(seq)
    original_seq = seq
    warnings: list[str] = []
    failure_trace: list[dict[str, Any]] = []
    seq_type_norm = (seq_type or "CDS").strip().upper()

    if seq_type_norm == "CDS" and remove_stop:
        seq2, trimmed = trim_stop_codon(seq)
        if trimmed:
            warnings.append("Stop codon removed")
        seq = seq2

    if len(seq) < max_len * 2:
        warnings.append("Sequence is short; primer design may be suboptimal")

    homo_f = normalize_seq(homo_f)
    homo_r = normalize_seq(homo_r)
    if not seq:
        raise ValueError("Input sequence is empty after normalization.")

    if min_len <= 0:
        raise ValueError("min_len must be positive.")
    if max_len < min_len:
        raise ValueError("max_len must be >= min_len.")

    profile_norm, tm_params = _resolve_tm_profile(tm_profile, Na, Mg, dNTPs, primer_conc)
    tm_calc_params = {"method": "Biopython_Tm_NN", **tm_params}

    is_non_cds = seq_type_norm in {"PROMOTER", "OTHER"}
    hp_filter_on = bool(non_cds_forbid_homopolymer and is_non_cds)

    pref_len = cds_len_pref if seq_type_norm == "CDS" else non_cds_len_pref
    pref_lo = int(pref_len[0])
    pref_hi = int(pref_len[1])
    if pref_lo > pref_hi:
        pref_lo, pref_hi = pref_hi, pref_lo
    if length_mode != "adaptive_soft":
        warnings.append(f"Unknown length_mode={length_mode}, fallback to adaptive_soft")

    selected_pair = None
    rescue_tier = "tier0"
    selected_reason = "best_tm"
    chosen_tm_min = float(tm_min)
    chosen_tm_max = float(tm_max)

    if seq_type_norm == "PROMOTER":
        fw_stats = _empty_stats()
        rv_stats = _empty_stats()
        fw_best_by_start: list[dict[str, Any]] = []
        rv_best_by_end: list[dict[str, Any]] = []
        n = len(seq)

        for i in range(0, n):
            cands, st = _iter_forward_candidates_at_start(
                seq=seq,
                start=i,
                min_len=min_len,
                max_len=max_len,
                tm_min=tm_min,
                tm_max=tm_max,
                target_tm=target_tm,
                gc_min=gc_min,
                gc_max=gc_max,
                pref_len_lo=pref_lo,
                pref_len_hi=pref_hi,
                tm_priority_weight=tm_priority_weight,
                forbid_homopolymer=hp_filter_on,
                tm_params=tm_params,
            )
            _merge_stats(fw_stats, st)
            best = _pick_best_by_score(cands)
            if best is not None:
                fw_best_by_start.append(best)

        for j in range(1, n + 1):
            cands, st = _iter_reverse_candidates_at_end(
                seq=seq,
                end=j,
                min_len=min_len,
                max_len=max_len,
                tm_min=tm_min,
                tm_max=tm_max,
                target_tm=target_tm,
                gc_min=gc_min,
                gc_max=gc_max,
                pref_len_lo=pref_lo,
                pref_len_hi=pref_hi,
                tm_priority_weight=tm_priority_weight,
                forbid_homopolymer=hp_filter_on,
                tm_params=tm_params,
            )
            _merge_stats(rv_stats, st)
            best = _pick_best_by_score(cands)
            if best is not None:
                rv_best_by_end.append(best)

        selected_pair = _select_pair(
            forwards=fw_best_by_start,
            reverses=rv_best_by_end,
            tm_pair_delta_max=float(tm_pair_delta_max),
            tm_delta_mode=tm_delta_mode,
            prefer_longest=True,
        )
        if selected_pair is None:
            reason = "; ".join(
                [
                    _build_failure_note("promoter forward", fw_stats, tm_min, tm_max, gc_min, gc_max),
                    _build_failure_note("promoter reverse", rv_stats, tm_min, tm_max, gc_min, gc_max),
                    f"tm_delta_mode={tm_delta_mode}, tm_pair_delta_max={tm_pair_delta_max}",
                ]
            )
            raise ValueError(reason)
        selected_reason = "best_tm_with_noncds_filter_longest_amplicon" if hp_filter_on else "best_tm_longest_amplicon"

    else:
        # CDS and Other keep end-anchored amplification, but choose at pair level.
        n = len(seq)
        attempts = [{"tier": "tier0", "tm_min": float(tm_min), "tm_max": float(tm_max), "max_len": int(max_len)}]
        if seq_type_norm == "CDS" and cds_recovery_enable:
            attempts = _build_cds_attempts(tm_min=float(tm_min), tm_max=float(tm_max), max_len=int(max_len))

        for a in attempts:
            a_tm_min = float(a["tm_min"])
            a_tm_max = float(a["tm_max"])
            a_max_len = int(a["max_len"])

            fws, fw_stats = _iter_forward_candidates_at_start(
                seq=seq,
                start=0,
                min_len=min_len,
                max_len=a_max_len,
                tm_min=a_tm_min,
                tm_max=a_tm_max,
                target_tm=target_tm,
                gc_min=gc_min,
                gc_max=gc_max,
                pref_len_lo=pref_lo,
                pref_len_hi=pref_hi,
                tm_priority_weight=tm_priority_weight,
                forbid_homopolymer=hp_filter_on,
                tm_params=tm_params,
            )
            rvs, rv_stats = _iter_reverse_candidates_at_end(
                seq=seq,
                end=n,
                min_len=min_len,
                max_len=a_max_len,
                tm_min=a_tm_min,
                tm_max=a_tm_max,
                target_tm=target_tm,
                gc_min=gc_min,
                gc_max=gc_max,
                pref_len_lo=pref_lo,
                pref_len_hi=pref_hi,
                tm_priority_weight=tm_priority_weight,
                forbid_homopolymer=hp_filter_on,
                tm_params=tm_params,
            )

            pair = _select_pair(
                forwards=fws,
                reverses=rvs,
                tm_pair_delta_max=float(tm_pair_delta_max),
                tm_delta_mode=tm_delta_mode,
                prefer_longest=False,
            )
            if pair is not None:
                selected_pair = pair
                rescue_tier = str(a["tier"])
                chosen_tm_min = a_tm_min
                chosen_tm_max = a_tm_max
                break

            failure_trace.append(
                {
                    "tier": str(a["tier"]),
                    "tm_min": a_tm_min,
                    "tm_max": a_tm_max,
                    "max_len": a_max_len,
                    "forward": dict(fw_stats),
                    "reverse": dict(rv_stats),
                    "reason": "No valid end-anchored pair after constraints.",
                }
            )

        if selected_pair is None:
            raise ValueError(
                "No valid primer pair found under constraints; "
                + "; ".join(
                    [
                        f"{x['tier']}[Tm {x['tm_min']:.1f}-{x['tm_max']:.1f}, max_len={x['max_len']}]"
                        f" Fpass={x['forward']['pass']} Rpass={x['reverse']['pass']}"
                        for x in failure_trace
                    ]
                )
            )

        selected_reason = "best_tm_with_noncds_filter" if hp_filter_on else "best_tm"

    fw = selected_pair["fw"]
    rv = selected_pair["rv"]
    tm_f = float(fw["tm"])
    tm_r = float(rv["tm"])
    tm_delta = float(selected_pair["tm_delta"])
    amplicon_len = int(selected_pair["amplicon_len"]) if seq_type_norm == "PROMOTER" else int(len(seq))
    gene_bind_f = str(fw["primer_seq"])
    gene_bind_r = str(rv["primer_seq"])
    bind_start_f = int(fw["bind_start"])
    bind_end_f = int(fw["bind_end"])
    bind_start_r = int(rv["bind_start"])
    bind_end_r = int(rv["bind_end"])
    amplicon_seq = seq[bind_start_f:bind_end_r] if seq_type_norm in {"PROMOTER", "OTHER"} else ""
    amplicon_len = int(bind_end_r - bind_start_f)
    gc_f = float(fw["gc"])
    gc_r = float(rv["gc"])

    if selected_pair.get("delta_relaxed"):
        warnings.append(f"Tm delta relaxed ({tm_delta:.2f} > {tm_pair_delta_max:.2f})")
        selected_reason = f"{selected_reason}|tm_delta_relaxed"
    if seq_type_norm == "PROMOTER":
        warnings.append(f"Promoter mode: selected longest amplicon {amplicon_len} bp")
    if rescue_tier != "tier0":
        warnings.append(f"CDS rescue applied: {rescue_tier} (Tm window {chosen_tm_min:.1f}-{chosen_tm_max:.1f})")
        selected_reason = f"{selected_reason}|{rescue_tier}"
    if tm_delta > float(tm_pair_delta_max):
        warnings.append(f"Tm delta warning ({tm_delta:.2f})")

    frame_check_warnings: list[str] = []
    if enable_frame_check:
        frame_check_warnings = _frame_check_messages(
            seq=seq,
            seq_type_norm=seq_type_norm,
            construct_type=construct_type,
            remove_stop=remove_stop,
            restriction_site_f=restriction_site_f,
            restriction_site_r=restriction_site_r,
            junction_extra_bases_f=junction_extra_bases_f,
            junction_extra_bases_r=junction_extra_bases_r,
        )
        warnings.extend(frame_check_warnings)

    final_f = homo_f + gene_bind_f
    final_r = homo_r + gene_bind_r

    return {
        "gene_name": gene_name,
        "gene_bind_f": gene_bind_f,
        "gene_bind_r": gene_bind_r,
        "bind_start_f": bind_start_f,
        "bind_end_f": bind_end_f,
        "bind_start_r": bind_start_r,
        "bind_end_r": bind_end_r,
        "amplicon_seq": amplicon_seq,
        "final_f": final_f,
        "final_r": final_r,
        "tm_f": tm_f,
        "tm_r": tm_r,
        "tm_delta": tm_delta,
        "gc_f": gc_f,
        "gc_r": gc_r,
        "amplicon_len": amplicon_len,
        "insert_seq": seq,
        "original_seq": original_seq,
        "warnings": "; ".join(warnings),
        "frame_check_warnings": "; ".join(frame_check_warnings),
        "seq_type": seq_type_norm,
        "construct_type": construct_type,
        "remove_stop": bool(remove_stop),
        "selected_reason": selected_reason,
        "rescue_tier": rescue_tier,
        "tm_profile": profile_norm,
        "tm_calc_params": tm_calc_params,
        "failure_trace": failure_trace,
    }
