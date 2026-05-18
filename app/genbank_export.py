import io
import re
import zipfile
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from primer_design import normalize_seq, revcomp


def _safe_name(s: str) -> str:
    s = (s or "").strip()
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_.-]", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "sample"


def _feature(start: int, end: int, ftype: str, **qualifiers: Any) -> SeqFeature:
    q = {k: [str(v)] for k, v in qualifiers.items() if v is not None and str(v) != ""}
    return SeqFeature(FeatureLocation(int(start), int(end), strand=1), type=ftype, qualifiers=q)


def _find_subseq(haystack: str, needle: str) -> tuple[int, int] | None:
    if not needle:
        return None
    i = haystack.find(needle)
    if i < 0:
        return None
    return i, i + len(needle)


def _make_construct(row: dict[str, Any]) -> tuple[str, dict[str, int], list[str]]:
    warnings: list[str] = []
    vector = normalize_seq(row.get("vector_sequence", "") or "")
    insert_seq = normalize_seq(row.get("insert_seq", "") or "")
    homo_f = normalize_seq(row.get("homo_f", "") or "")
    homo_r = normalize_seq(row.get("homo_r", "") or "")
    rsf = normalize_seq(row.get("restriction_site_f", "") or "")
    rsr = normalize_seq(row.get("restriction_site_r", "") or "")
    jef = max(0, int(row.get("junction_extra_bases_f") or 0))
    jer = max(0, int(row.get("junction_extra_bases_r") or 0))
    pad_f = "N" * jef
    pad_r = "N" * jer

    if vector and homo_f and homo_r:
        start_f = vector.find(homo_f)
        start_r = vector.rfind(homo_r)
        if start_f >= 0 and start_r >= 0:
            end_f = start_f + len(homo_f)
            if end_f <= start_r:
                seq = vector[:end_f] + insert_seq + vector[start_r:]
                return seq, {"insert_start": end_f, "insert_end": end_f + len(insert_seq)}, warnings
        warnings.append("Warning: HOMO_F/HOMO_R not locatable in vector_sequence; fallback to minimal construct.")

    seq = vector + rsf + pad_f + insert_seq + pad_r + rsr
    insert_start = len(vector) + len(rsf) + len(pad_f)
    insert_end = insert_start + len(insert_seq)
    if not vector:
        warnings.append("Warning: vector_sequence missing; exported minimal construct only.")
    return seq, {"insert_start": insert_start, "insert_end": insert_end}, warnings


def _build_record(row: dict[str, Any]) -> SeqRecord:
    seq, pos, extra_warnings = _make_construct(row)
    backbone_name = str(row.get("backbone_name", "")).strip() or "backbone"
    gene_name = str(row.get("gene_name", "")).strip() or "gene"
    seq_type = str(row.get("seq_type", "")).strip() or "CDS"
    construct_id = f"{_safe_name(backbone_name)}_{_safe_name(gene_name)}_construct"
    record = SeqRecord(Seq(seq), id=construct_id, name=construct_id[:16], description=f"{backbone_name} + {gene_name}")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "linear"

    insert_start = int(pos["insert_start"])
    insert_end = int(pos["insert_end"])
    n = len(seq)
    rsf = normalize_seq(row.get("restriction_site_f", "") or "")
    rsr = normalize_seq(row.get("restriction_site_r", "") or "")
    jef = max(0, int(row.get("junction_extra_bases_f") or 0))
    jer = max(0, int(row.get("junction_extra_bases_r") or 0))

    record.features.append(_feature(0, max(1, n), "source", organism="synthetic construct"))
    record.features.append(_feature(0, max(1, n), "misc_feature", label="vector/backbone", note=backbone_name))
    if insert_end > insert_start:
        gene_type = "CDS" if seq_type.upper() == "CDS" else "misc_feature"
        record.features.append(_feature(insert_start, insert_end, gene_type, label="inserted gene/CDS", gene=gene_name))

    homo_f = normalize_seq(row.get("homo_f", "") or "")
    homo_r = normalize_seq(row.get("homo_r", "") or "")
    if homo_f:
        loc = _find_subseq(seq, homo_f)
        if loc:
            record.features.append(_feature(loc[0], loc[1], "misc_feature", label="HOMO_F", sequence=homo_f))
    if homo_r:
        loc = _find_subseq(seq, homo_r)
        if loc:
            record.features.append(_feature(loc[0], loc[1], "misc_feature", label="HOMO_R", sequence=homo_r))

    final_f = normalize_seq(row.get("final_f", "") or "")
    final_r = normalize_seq(row.get("final_r", "") or "")
    if final_f:
        loc = _find_subseq(seq, final_f)
        if loc is None:
            loc = (insert_start, min(insert_start + max(1, len(final_f)), n))
        record.features.append(_feature(loc[0], loc[1], "primer_bind", label="final_F_primer", sequence=final_f))
    if final_r:
        loc = _find_subseq(seq, revcomp(final_r))
        if loc is None:
            loc = (max(0, insert_end - max(1, len(final_r))), min(insert_end, n))
        record.features.append(_feature(loc[0], loc[1], "primer_bind", label="final_R_primer", sequence=final_r))

    rsf_start = insert_start - (len(rsf) + jef)
    if rsf:
        record.features.append(
            _feature(max(0, rsf_start), max(0, rsf_start + len(rsf)), "misc_feature", label="restriction_site_F", sequence=rsf)
        )
    rsr_start = insert_end + jer
    if rsr:
        record.features.append(
            _feature(max(0, rsr_start), max(0, rsr_start + len(rsr)), "misc_feature", label="restriction_site_R", sequence=rsr)
        )

    j5_start = max(0, insert_start - 1)
    j5_end = min(n, insert_start + 1)
    j3_start = max(0, insert_end - 1)
    j3_end = min(n, insert_end + 1)
    record.features.append(_feature(j5_start, max(j5_start + 1, j5_end), "misc_feature", label="junction_5prime"))
    record.features.append(_feature(j3_start, max(j3_start + 1, j3_end), "misc_feature", label="junction_3prime"))

    note = str(row.get("vector_map_note", "") or "").strip()
    warnings = str(row.get("warnings", "") or "").strip()
    merged = "; ".join([x for x in [warnings, *extra_warnings, note] if x])
    if merged:
        record.features.append(_feature(0, max(1, n), "misc_feature", label="design_notes", note=merged))
    return record


def _record_to_gbk_bytes(record: SeqRecord) -> bytes:
    sio = io.StringIO()
    SeqIO.write(record, sio, "genbank")
    return sio.getvalue().encode("utf-8")


def build_genbank_bundle(rows: list[dict[str, Any]]) -> tuple[str, bytes, str]:
    if not rows:
        raise ValueError("No rows to export as GenBank.")
    files: list[tuple[str, bytes]] = []
    for row in rows:
        record = _build_record(row)
        file_name = f"{_safe_name(row.get('backbone_name', 'backbone'))}_{_safe_name(row.get('gene_name', 'gene'))}_construct.gbk"
        files.append((file_name, _record_to_gbk_bytes(record)))

    if len(files) == 1:
        return files[0][0], files[0][1], "application/octet-stream"

    bio = io.BytesIO()
    with zipfile.ZipFile(bio, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, data in files:
            zf.writestr(name, data)
    return "Construct_Maps.zip", bio.getvalue(), "application/zip"

