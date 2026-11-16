from typing import Any, List, Tuple
import io
import re
from Bio import SeqIO
from primer_design import normalize_seq


def sanitize_name(raw: str) -> str:
    s = (raw or "").strip()
    s = s.replace("\ufeff", "")
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_.-]", "_", s)
    s = re.sub(r"_+", "_", s).strip("_.")
    return s or "Gene"


def derive_gene_name(raw: str) -> str:
    """Derive a short, user-friendly gene name like OsTF1 / ZmTF1 / TF1."""
    s = sanitize_name(raw)
    s = re.sub(r"\.(fasta|fa|fna|txt|gb|gbk|genbank|dna)$", "", s, flags=re.IGNORECASE)
    token = s.split()[0].strip()

    m = re.fullmatch(r"(?i)([a-z]{0,4})tf0*([0-9]+)", token)
    if m:
        prefix = m.group(1) or ""
        num = str(int(m.group(2)))
        if prefix:
            prefix = prefix[0].upper() + prefix[1:].lower()
        return f"{prefix}TF{num}"
    return sanitize_name(token)


def parse_fasta_text(text: str) -> List[Tuple[str, str]]:
    text = (text or "").strip()
    if not text:
        return []
    if text.startswith(">"):
        records = []
        name = None
        seq_lines: List[str] = []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((derive_gene_name(name), normalize_seq("".join(seq_lines))))
                name = line[1:].strip().split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if name is not None:
            records.append((derive_gene_name(name), normalize_seq("".join(seq_lines))))
        return records
    return [("Gene1", normalize_seq(text))]


def file_base_name(filename: str) -> str:
    return sanitize_name(re.sub(r"\.[^.]+$", "", filename or ""))


def _is_snapgene_magic(content: bytes) -> bool:
    if len(content) < 14:
        return False
    return content[:5] == b"\t\x00\x00\x00\x0e" and content[5:13] == b"SnapGene"


def _looks_binary(content: bytes) -> bool:
    if not content:
        return False
    sample = content[:1024]
    if b"\x00" in sample:
        return True
    non_text = sum(1 for b in sample if b < 9 or (b > 13 and b < 32))
    return (non_text / max(1, len(sample))) > 0.10


def _feature_to_json(feature) -> dict[str, Any]:
    loc = feature.location
    parts = getattr(loc, "parts", None) or [loc]
    out_parts = []
    for part in parts:
        out_parts.append(
            {
                "start": int(part.start),
                "end": int(part.end),
                "strand": int(part.strand or 0),
            }
        )
    return {
        "type": str(feature.type or "misc_feature"),
        "parts": out_parts,
        "qualifiers": {
            str(k): [str(x) for x in (v if isinstance(v, list) else [v])]
            for k, v in (feature.qualifiers or {}).items()
        },
    }


def _snapgene_reader_feature_to_json(feature: dict[str, Any]) -> dict[str, Any] | None:
    try:
        start = int(feature.get("start"))
        end = int(feature.get("end"))
    except Exception:
        return None
    strand_raw = str(feature.get("strand") or ".")
    strand = -1 if strand_raw == "-" else (1 if strand_raw == "+" else 0)
    qualifiers = feature.get("qualifiers") or {}
    fixed_qualifiers: dict[str, list[str]] = {}
    for key, value in qualifiers.items():
        if isinstance(value, list):
            fixed_qualifiers[str(key)] = [str(x) for x in value]
        else:
            fixed_qualifiers[str(key)] = [str(value)]
    label = feature.get("name")
    if label and "label" not in fixed_qualifiers:
        fixed_qualifiers["label"] = [str(label)]
    return {
        "type": str(feature.get("type") or "misc_feature"),
        "parts": [{"start": start, "end": end, "strand": strand}],
        "qualifiers": fixed_qualifiers,
    }


def _record_to_parsed(rec, fallback_name: str, source_format: str) -> dict[str, Any] | None:
    seq = normalize_seq(str(rec.seq))
    if not seq:
        return None
    rec_name = rec.name or rec.id or fallback_name
    topology = str(rec.annotations.get("topology") or "").lower()
    if topology not in {"linear", "circular"}:
        topology = ""
    return {
        "name": derive_gene_name(rec_name),
        "sequence": seq,
        "features": [_feature_to_json(f) for f in rec.features],
        "topology": topology,
        "source_format": source_format,
        "record_name": str(rec_name),
    }


def _parse_dna_content_detailed(content: bytes, filename: str) -> list[dict[str, Any]]:
    parse_errors = []

    # Parser 1: Biopython preserves SnapGene features and topology well.
    try:
        bio = io.BytesIO(content)
        records = []
        for rec in SeqIO.parse(bio, "snapgene"):
            parsed = _record_to_parsed(rec, file_base_name(filename), "snapgene")
            if parsed:
                records.append(parsed)
        if records:
            return records
        parse_errors.append("Bio.SeqIO(snapgene): no readable records")
    except Exception as e:
        parse_errors.append(f"Bio.SeqIO(snapgene) parse failed: {e}")

    # Parser 2: snapgene-reader sequence/features fallback.
    try:
        from snapgene_reader import snapgene_file_to_dict  # type: ignore

        bio = io.BytesIO(content)
        d = snapgene_file_to_dict(fileobject=bio)
        seq = normalize_seq(d.get("seq", "") or "")
        if seq:
            rec_name = d.get("name") or file_base_name(filename)
            features = []
            for item in d.get("features") or []:
                parsed_feature = _snapgene_reader_feature_to_json(item)
                if parsed_feature:
                    features.append(parsed_feature)
            return [
                {
                    "name": derive_gene_name(rec_name),
                    "sequence": seq,
                    "features": features,
                    "topology": "",
                    "source_format": "snapgene",
                    "record_name": str(rec_name),
                }
            ]
        parse_errors.append("snapgene-reader: no sequence field in file")
    except ImportError:
        parse_errors.append("snapgene-reader dependency missing")
    except Exception as e:
        parse_errors.append(f"snapgene-reader parse failed: {e}")

    raise ValueError(
        "Failed to parse SnapGene DNA content. "
        "Possible causes: dependency missing, unsupported SnapGene variant, or file corruption. "
        f"Details: {' | '.join(parse_errors)}"
    )


def parse_uploaded_records(file) -> list[dict[str, Any]]:
    """
    Supports:
      - FASTA/TXT (single or multi-FASTA)
      - GenBank (.gb/.gbk/.genbank)
      - SnapGene (.dna) when snapgene-reader is available
    """
    name = file.name or "uploaded_sequence"
    lower = name.lower()
    content = file.read()

    # Content sniffing: some SnapGene files are mis-labeled as .fa/.txt.
    if _is_snapgene_magic(content):
        return _parse_dna_content_detailed(content, name)

    if lower.endswith((".fasta", ".fa", ".fna", ".txt")):
        if _looks_binary(content):
            raise ValueError(
                "This file looks binary, not FASTA text. "
                "If it is a SnapGene file, keep extension as .dna or export FASTA from SnapGene."
            )
        text = content.decode("utf-8", errors="ignore").strip()
        if not text:
            raise ValueError("Empty text file.")
        if text.startswith(">"):
            recs = parse_fasta_text(text)
            if not recs:
                raise ValueError("No FASTA records found.")
            return [
                {
                    "name": rec_name,
                    "sequence": seq,
                    "features": [],
                    "topology": "",
                    "source_format": "fasta",
                    "record_name": rec_name,
                }
                for rec_name, seq in recs
            ]
        raw = re.sub(r"\s+", "", text).upper()
        if raw and re.fullmatch(r"[ACGTRYSWKMBDHVN]+", raw):
            rec_name = derive_gene_name(file_base_name(name))
            return [
                {
                    "name": rec_name,
                    "sequence": normalize_seq(raw),
                    "features": [],
                    "topology": "",
                    "source_format": "text",
                    "record_name": rec_name,
                }
            ]
        raise ValueError("Text file is not valid FASTA or plain DNA sequence.")

    if lower.endswith((".gb", ".gbk", ".genbank")):
        bio = io.StringIO(content.decode("utf-8-sig"))
        records = []
        for rec in SeqIO.parse(bio, "genbank"):
            parsed = _record_to_parsed(rec, file_base_name(name), "genbank")
            if parsed:
                records.append(parsed)
        return records

    if lower.endswith(".dna"):
        return _parse_dna_content_detailed(content, name)

    raise ValueError(f"Unsupported file type: {name}")


def parse_uploaded_file(file) -> List[Tuple[str, str]]:
    return [(r["name"], r["sequence"]) for r in parse_uploaded_records(file)]
