from typing import List, Tuple
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


def _parse_dna_content(content: bytes, filename: str) -> List[Tuple[str, str]]:
    parse_errors = []

    # Parser 1: snapgene-reader
    try:
        from snapgene_reader import snapgene_file_to_dict  # type: ignore

        bio = io.BytesIO(content)
        d = snapgene_file_to_dict(fileobject=bio)
        seq = normalize_seq(d.get("seq", "") or "")
        if seq:
            rec_name = d.get("name") or file_base_name(filename)
            return [(derive_gene_name(rec_name), seq)]
        parse_errors.append("snapgene-reader: no sequence field in file")
    except ImportError:
        parse_errors.append("snapgene-reader dependency missing")
    except Exception as e:
        parse_errors.append(f"snapgene-reader parse failed: {e}")

    # Parser 2: Biopython snapgene parser fallback
    try:
        bio = io.BytesIO(content)
        records = []
        for rec in SeqIO.parse(bio, "snapgene"):
            rec_name = rec.name or rec.id or file_base_name(filename)
            seq = normalize_seq(str(rec.seq))
            if seq:
                records.append((derive_gene_name(rec_name), seq))
        if records:
            return records
        parse_errors.append("Bio.SeqIO(snapgene): no readable records")
    except Exception as e:
        parse_errors.append(f"Bio.SeqIO(snapgene) parse failed: {e}")

    raise ValueError(
        "Failed to parse SnapGene DNA content. "
        "Possible causes: dependency missing, unsupported SnapGene variant, or file corruption. "
        f"Details: {' | '.join(parse_errors)}"
    )


def parse_uploaded_file(file) -> List[Tuple[str, str]]:
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
        return _parse_dna_content(content, name)

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
            return recs
        raw = re.sub(r"\s+", "", text).upper()
        if raw and re.fullmatch(r"[ACGTRYSWKMBDHVN]+", raw):
            return [(derive_gene_name(file_base_name(name)), normalize_seq(raw))]
        raise ValueError("Text file is not valid FASTA or plain DNA sequence.")

    if lower.endswith((".gb", ".gbk", ".genbank")):
        bio = io.BytesIO(content)
        records = []
        for rec in SeqIO.parse(bio, "genbank"):
            rec_name = rec.name or rec.id or file_base_name(name)
            records.append((derive_gene_name(rec_name), normalize_seq(str(rec.seq))))
        return records

    if lower.endswith(".dna"):
        return _parse_dna_content(content, name)

    raise ValueError(f"Unsupported file type: {name}")
