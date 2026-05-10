import time


def _yn(v) -> str:
    if isinstance(v, str):
        return v
    return "Yes" if bool(v) else "No"


def _fmt_float(v) -> str:
    if v is None or str(v) == "":
        return "NA"
    try:
        return f"{float(v):.2f}"
    except Exception:
        return "NA"


def build_fragment_length_txt(rows) -> bytes:
    lines = []
    lines.append("# Fragment Length Record")
    lines.append(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append(
        "Backbone\tGene\tSeqType\tRemoveStop\tAmplicon_bp\tTmF\tTmR\tTmDelta\tTmProfile\tSelectedReason\tPrimerF\tPrimerR\tCorePrimerF_noHomology\tCorePrimerR_noHomology"
    )
    for r in rows:
        amp = r.get("amplicon_len")
        amp_s = str(int(amp)) if amp is not None and str(amp) != "" else "NA"
        tm_f_s = _fmt_float(r.get("tm_f"))
        tm_r_s = _fmt_float(r.get("tm_r"))
        tm_d_s = _fmt_float(r.get("tm_delta"))
        lines.append(
            "\t".join(
                [
                    str(r.get("backbone_name", "")),
                    str(r.get("gene_name", "")),
                    str(r.get("seq_type", "")),
                    _yn(r.get("remove_stop", "")),
                    amp_s,
                    tm_f_s,
                    tm_r_s,
                    tm_d_s,
                    str(r.get("tm_profile", "")),
                    str(r.get("selected_reason", "")),
                    str(r.get("primer_name_f", "")),
                    str(r.get("primer_name_r", "")),
                    str(r.get("gene_bind_f", "")),
                    str(r.get("gene_bind_r", "")),
                ]
            )
        )
    text = "\n".join(lines) + "\n"
    # UTF-8 with BOM makes it easier to open directly in Windows Notepad/Excel.
    return text.encode("utf-8-sig")


def build_fragment_length_simple_txt(rows) -> bytes:
    lines = []
    lines.append("# Fragment Length (Simple)")
    lines.append("Sample\tAmplicon_bp\tTmF\tTmR\tTmDelta\tCorePrimerF\tCorePrimerR")
    for r in rows:
        amp = r.get("amplicon_len")
        amp_s = str(int(amp)) if amp is not None and str(amp) != "" else "NA"
        tm_f_s = _fmt_float(r.get("tm_f"))
        tm_r_s = _fmt_float(r.get("tm_r"))
        tm_d_s = _fmt_float(r.get("tm_delta"))
        backbone = str(r.get("backbone_name", "")).strip()
        gene = str(r.get("gene_name", "")).strip()
        sample = "_".join([x for x in [backbone, gene] if x]) or str(r.get("primer_name_f", "")).replace("-F", "")
        lines.append(f"{sample}\t{amp_s}\t{tm_f_s}\t{tm_r_s}\t{tm_d_s}\t{str(r.get('gene_bind_f', ''))}\t{str(r.get('gene_bind_r', ''))}")
    text = "\n".join(lines) + "\n"
    return text.encode("utf-8-sig")


def build_non_cds_amplicon_fasta(rows) -> bytes:
    lines = []
    count = 0
    for r in rows:
        seq_type = str(r.get("seq_type", "")).strip().upper()
        if seq_type not in {"PROMOTER", "OTHER"}:
            continue
        amp_seq = str(r.get("amplicon_seq", "")).strip().upper()
        if not amp_seq:
            continue
        count += 1
        backbone = str(r.get("backbone_name", "")).strip()
        gene = str(r.get("gene_name", "")).strip()
        sample = "_".join([x for x in [backbone, gene] if x]) or f"non_cds_{count}"
        header = (
            f">{sample}"
            f"|seq_type={seq_type}"
            f"|amp_bp={str(r.get('amplicon_len', 'NA'))}"
            f"|tmF={_fmt_float(r.get('tm_f'))}"
            f"|tmR={_fmt_float(r.get('tm_r'))}"
            f"|dTm={_fmt_float(r.get('tm_delta'))}"
        )
        lines.append(header)
        lines.append(amp_seq)
    if count == 0:
        return b""
    text = "\n".join(lines) + "\n"
    return text.encode("utf-8-sig")
