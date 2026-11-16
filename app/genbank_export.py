import io
import json
import re
import zipfile
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from primer_design import normalize_seq, revcomp


KNOWN_ENZYMES = {
    "GGATCC": "BamHI",
    "GAATTC": "EcoRI",
    "AAGCTT": "HindIII",
    "TCTAGA": "XbaI",
    "ACTAGT": "SpeI",
    "GGTACC": "KpnI",
    "GTCGAC": "SalI",
    "CCCGGG": "SmaI",
}


def _safe_name(s: str) -> str:
    s = (s or "").strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_.-]", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "sample"


def _qualifier_values(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(x) for x in value if str(x) != ""]
    return [str(value)] if str(value) != "" else []


def _feature(start: int, end: int, ftype: str, strand: int | None = 1, **qualifiers: Any) -> SeqFeature:
    start = max(0, int(start))
    end = max(start + 1, int(end))
    q = {k: _qualifier_values(v) for k, v in qualifiers.items()}
    q = {k: v for k, v in q.items() if v}
    return SeqFeature(FeatureLocation(start, end, strand=strand), type=ftype, qualifiers=q)


def _find_all(seq: str, needle: str) -> list[int]:
    seq = normalize_seq(seq)
    needle = normalize_seq(needle)
    if not seq or not needle:
        return []
    return [i for i in range(0, len(seq) - len(needle) + 1) if seq.startswith(needle, i)]


def _unique_position(seq: str, needle: str, label: str) -> int:
    positions = _find_all(seq, needle)
    if len(positions) == 1:
        return positions[0]
    if len(positions) == 0:
        raise ValueError(f"{label} 位点在骨架中找到 0 个，请检查骨架序列或切点设置。")
    raise ValueError(f"{label} 位点找到多个（{len(positions)} 个），请进入高级模式手动选择切点。")


def _enzyme_label(seq: str, fallback: str) -> str:
    seq = normalize_seq(seq)
    return KNOWN_ENZYMES.get(seq, fallback)


def _enzyme_labels(row: dict[str, Any], rsf: str, rsr: str) -> tuple[str, str]:
    labels = [x.strip() for x in str(row.get("enzyme_sites") or "").split("/") if x.strip()]
    f_label = labels[0] if len(labels) >= 1 else _enzyme_label(rsf, "restriction_site_F")
    r_label = labels[1] if len(labels) >= 2 else _enzyme_label(rsr, "restriction_site_R")
    return f_label, r_label


def _resolve_insert_feature_type(row: dict[str, Any]) -> str:
    explicit = str(row.get("insert_feature_type") or row.get("map_insert_feature_type") or "").strip()
    if explicit in {"CDS", "promoter", "misc_feature"}:
        return explicit
    seq_type = str(row.get("seq_type", "")).strip().lower()
    if seq_type == "cds":
        return "CDS"
    if seq_type == "promoter":
        return "promoter"
    return "misc_feature"


def _resolve_homology_feature_type(row: dict[str, Any]) -> str:
    explicit = str(row.get("homology_feature_type") or row.get("map_homology_feature_type") or "").strip()
    if explicit in {"misc_feature", "primer_bind"}:
        return explicit
    return "misc_feature"


def _resolve_insert_name(row: dict[str, Any]) -> str:
    return str(row.get("insert_feature_name") or row.get("map_insert_feature_name") or row.get("gene_name") or "insert").strip()


def _resolve_cut_window(row: dict[str, Any], vector: str) -> tuple[int, int, list[str]]:
    warnings: list[str] = []
    homo_f = normalize_seq(row.get("homo_f", "") or "")
    homo_r = normalize_seq(row.get("homo_r", "") or "")
    if homo_f and homo_r:
        f_positions = _find_all(vector, homo_f)
        r_seq = revcomp(homo_r)
        r_positions = _find_all(vector, r_seq)
        if len(f_positions) == 1 and len(r_positions) == 1:
            replace_start = f_positions[0] + len(homo_f)
            replace_end = r_positions[0]
            if replace_start <= replace_end:
                return replace_start, replace_end, warnings
            raise ValueError("无法确认插入方向，请检查同源臂是否与插入片段末端匹配。")
        warnings.append(
            f"同源臂未能唯一定位：F={len(f_positions)} 个，R={len(r_positions)} 个；改用酶切位点定位。"
        )

    rsf = normalize_seq(row.get("restriction_site_f", "") or "")
    rsr = normalize_seq(row.get("restriction_site_r", "") or "")
    if rsf and rsr:
        f_label, r_label = _enzyme_labels(row, rsf, rsr)
        f_pos = _unique_position(vector, rsf, f_label)
        r_pos = _unique_position(vector, rsr, r_label)
        replace_start = f_pos + len(rsf)
        replace_end = r_pos
        if replace_start <= replace_end:
            return replace_start, replace_end, warnings
        raise ValueError("无法确认插入方向，请检查同源臂是否与插入片段末端匹配。")

    raise ValueError("无法自动定位插入切点，请在高级模式提供同源臂或限制性酶切位点。")


def _map_part(start: int, end: int, replace_start: int, replace_end: int, insert_len: int) -> tuple[int, int] | None:
    removed_len = replace_end - replace_start
    delta = insert_len - removed_len
    if end <= replace_start:
        return start, end
    if start >= replace_end:
        return start + delta, end + delta
    if start < replace_start and end > replace_end:
        return start, end + delta
    if start < replace_start < end <= replace_end:
        return start, replace_start
    if replace_start <= start < replace_end < end:
        return replace_start + insert_len, end + delta
    return None


def _load_vector_features(row: dict[str, Any]) -> list[dict[str, Any]]:
    raw = row.get("vector_features_json") or row.get("backbone_features_json") or "[]"
    if isinstance(raw, list):
        return raw
    try:
        data = json.loads(str(raw or "[]"))
        return data if isinstance(data, list) else []
    except Exception:
        return []


def _copy_qualifiers(item: dict[str, Any]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    qualifiers = item.get("qualifiers") or {}
    if isinstance(qualifiers, dict):
        for key, value in qualifiers.items():
            vals = _qualifier_values(value)
            if vals:
                out[str(key)] = vals
    return out


def _remap_vector_features(
    features: list[dict[str, Any]],
    replace_start: int,
    replace_end: int,
    insert_len: int,
) -> tuple[list[SeqFeature], int]:
    out: list[SeqFeature] = []
    dropped = 0
    for item in features:
        if str(item.get("type") or "") == "source":
            continue
        parts = item.get("parts") or []
        mapped_parts = []
        for part in parts:
            try:
                mapped = _map_part(
                    int(part["start"]),
                    int(part["end"]),
                    replace_start,
                    replace_end,
                    insert_len,
                )
            except Exception:
                mapped = None
            if mapped and mapped[1] > mapped[0]:
                mapped_parts.append(FeatureLocation(mapped[0], mapped[1], strand=int(part.get("strand") or 0) or None))
        if not mapped_parts:
            dropped += 1
            continue
        location = mapped_parts[0] if len(mapped_parts) == 1 else CompoundLocation(mapped_parts)
        out.append(SeqFeature(location, type=str(item.get("type") or "misc_feature"), qualifiers=_copy_qualifiers(item)))
    return out, dropped


def _find_feature_location(seq: str, needle: str, label: str) -> tuple[int, int] | None:
    positions = _find_all(seq, needle)
    if len(positions) == 1:
        return positions[0], positions[0] + len(normalize_seq(needle))
    return None


def _make_construct(row: dict[str, Any]) -> tuple[str, dict[str, int], list[str], int, int]:
    warnings: list[str] = []
    vector = normalize_seq(row.get("vector_sequence", "") or "")
    insert_seq = normalize_seq(row.get("insert_seq", "") or "")
    if not insert_seq:
        raise ValueError("没有可导出的插入片段序列，请先成功生成引物。")
    if not vector:
        raise ValueError("骨架缺少完整序列，无法生成完整质粒图谱。")

    replace_start, replace_end, locate_warnings = _resolve_cut_window(row, vector)
    warnings.extend(locate_warnings)
    seq = vector[:replace_start] + insert_seq + vector[replace_end:]
    pos = {"insert_start": replace_start, "insert_end": replace_start + len(insert_seq)}
    return seq, pos, warnings, replace_start, replace_end


def build_genbank_record(row: dict[str, Any]) -> SeqRecord:
    seq, pos, extra_warnings, replace_start, replace_end = _make_construct(row)
    backbone_name = str(row.get("backbone_name", "")).strip() or "backbone"
    gene_name = str(row.get("gene_name", "")).strip() or "gene"
    insert_name = _resolve_insert_name(row)
    plasmid_name = str(row.get("plasmid_name") or "").strip()
    construct_id = _safe_name(plasmid_name) if plasmid_name else f"{_safe_name(backbone_name)}_{_safe_name(gene_name)}_construct"
    description = plasmid_name or f"{backbone_name} + {gene_name}"
    record = SeqRecord(Seq(seq), id=construct_id, name=construct_id[:16], description=description)
    record.annotations["molecule_type"] = "DNA"
    topology = str(row.get("vector_topology") or "").strip().lower()
    record.annotations["topology"] = topology if topology in {"linear", "circular"} else "circular"

    insert_start = int(pos["insert_start"])
    insert_end = int(pos["insert_end"])
    n = len(seq)
    insert_len = insert_end - insert_start

    record.features.append(_feature(0, max(1, n), "source", organism="synthetic construct", note=description))
    preserved, dropped_count = _remap_vector_features(_load_vector_features(row), replace_start, replace_end, insert_len)
    record.features.extend(preserved)

    insert_type = _resolve_insert_feature_type(row)
    insert_note = str(row.get("insert_feature_note") or row.get("map_insert_feature_note") or "Inserted sequence").strip()
    insert_qualifiers: dict[str, Any] = {"label": insert_name, "note": insert_note}
    if insert_type == "CDS":
        insert_qualifiers["gene"] = str(row.get("insert_gene") or insert_name)
        if row.get("insert_locus_tag"):
            insert_qualifiers["locus_tag"] = str(row["insert_locus_tag"])
        if row.get("insert_translation"):
            insert_qualifiers["translation"] = str(row["insert_translation"])
        insert_qualifiers["codon_start"] = "1"
    record.features.append(_feature(insert_start, insert_end, insert_type, **insert_qualifiers))

    homology_type = _resolve_homology_feature_type(row)
    homo_f = normalize_seq(row.get("homo_f", "") or "")
    homo_r = normalize_seq(row.get("homo_r", "") or "")
    if homo_f:
        loc = _find_feature_location(seq, homo_f, "F_homology_arm")
        if loc:
            record.features.append(_feature(loc[0], loc[1], homology_type, label="F_homology_arm", note="Forward homology arm"))
    if homo_r:
        r_seq = revcomp(homo_r)
        loc = _find_feature_location(seq, r_seq, "R_homology_arm")
        if loc:
            record.features.append(_feature(loc[0], loc[1], homology_type, strand=-1, label="R_homology_arm", note="Reverse homology arm"))

    rsf = normalize_seq(row.get("restriction_site_f", "") or "")
    rsr = normalize_seq(row.get("restriction_site_r", "") or "")
    f_label, r_label = _enzyme_labels(row, rsf, rsr)
    if rsf:
        loc = _find_feature_location(seq, rsf, f_label)
        fallback = (insert_start - len(rsf), insert_start)
        if (loc is None) and fallback[0] >= 0 and seq[fallback[0] : fallback[1]] == rsf:
            loc = fallback
        if loc:
            record.features.append(_feature(loc[0], loc[1], "misc_feature", label=f_label, note=f"{f_label} restriction site; direction=forward"))
    if rsr:
        loc = _find_feature_location(seq, rsr, r_label)
        fallback = (insert_end, insert_end + len(rsr))
        if (loc is None) and fallback[1] <= len(seq) and seq[fallback[0] : fallback[1]] == rsr:
            loc = fallback
        if loc:
            record.features.append(_feature(loc[0], loc[1], "misc_feature", strand=-1, label=r_label, note=f"{r_label} restriction site; direction=reverse"))

    final_f = normalize_seq(row.get("final_f", "") or "")
    final_r = normalize_seq(row.get("final_r", "") or "")
    if final_f:
        loc = _find_feature_location(seq, final_f, "final_F_primer")
        if loc:
            record.features.append(_feature(loc[0], loc[1], "primer_bind", label="final_F_primer", note="Final forward primer"))
    if final_r:
        loc = _find_feature_location(seq, revcomp(final_r), "final_R_primer")
        if loc:
            record.features.append(_feature(loc[0], loc[1], "primer_bind", strand=-1, label="final_R_primer", note="Final reverse primer"))

    notes = [str(row.get("warnings") or "").strip(), *extra_warnings]
    if dropped_count:
        notes.append(f"{dropped_count} original backbone feature(s) inside the replaced cloning interval were not retained.")
    note = str(row.get("vector_map_note", "") or "").strip()
    if note:
        notes.append(note)
    merged = "; ".join(x for x in notes if x)
    if merged:
        record.features.append(_feature(0, max(1, n), "misc_feature", label="design_notes", note=merged))
    return record


def _record_to_gbk_bytes(record: SeqRecord) -> bytes:
    sio = io.StringIO()
    SeqIO.write(record, sio, "genbank")
    return sio.getvalue().encode("utf-8")


def build_genbank_files(rows: list[dict[str, Any]]) -> list[tuple[str, bytes]]:
    if not rows:
        raise ValueError("No rows to export as GenBank.")
    files: list[tuple[str, bytes]] = []
    for row in rows:
        record = build_genbank_record(row)
        output_filename = str(row.get("output_filename") or "").strip()
        plasmid_name = str(row.get("plasmid_name") or "").strip()
        if output_filename:
            file_name = _safe_name(output_filename)
            if not file_name.lower().endswith(".gbk"):
                file_name += ".gbk"
        elif plasmid_name:
            file_name = f"{_safe_name(plasmid_name)}.gbk"
        else:
            file_name = f"{_safe_name(row.get('backbone_name', 'backbone'))}_{_safe_name(row.get('gene_name', 'gene'))}_construct.gbk"
        files.append((file_name, _record_to_gbk_bytes(record)))
    return files


def build_genbank_bundle(rows: list[dict[str, Any]]) -> tuple[str, bytes, str]:
    files = build_genbank_files(rows)
    if len(files) == 1:
        return files[0][0], files[0][1], "application/octet-stream"

    bio = io.BytesIO()
    with zipfile.ZipFile(bio, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, data in files:
            zf.writestr(name, data)
    return "Construct_Maps.zip", bio.getvalue(), "application/zip"
