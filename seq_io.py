from typing import List, Tuple
import re
from Bio import SeqIO
from primer_design import normalize_seq

def derive_gene_name(raw: str) -> str:
    """Derive a short, user-friendly gene name like OsTF1 / ZmTF1 / TF1.

    Rules:
      - Use the token before whitespace and before '-' (if present)
      - Match optional 0–4 letter prefix + 'TF' + number (case-insensitive)
        Examples:
          'TF1-LOC_Os09g...' -> 'TF1'
          'OsTF12-anything' -> 'OsTF12'
          'zmtf003' -> 'ZmTF3'
      - If no TF-style pattern is found, fall back to the base token.
    """
    s = (raw or "").strip()
    s = re.sub(r"\s+", " ", s)
    if not s:
        return s

    # remove common extensions if present
    s = re.sub(r"\.(fasta|fa|fna|txt|gb|gbk|genbank|dna)$", "", s, flags=re.IGNORECASE)

    token = s.split("-", 1)[0].strip()
    token = token.split()[0].strip()

    m = re.fullmatch(r"(?i)([a-z]{0,4})tf0*([0-9]+)", token)
    if m:
        prefix = m.group(1) or ""
        num = str(int(m.group(2)))
        if prefix:
            prefix = prefix[0].upper() + prefix[1:].lower()
        return f"{prefix}TF{num}"

    return token

def parse_fasta_text(text: str) -> List[Tuple[str,str]]:
    text = (text or "").strip()
    if not text:
        return []
    if text.startswith(">"):
        records = []
        name = None
        seq_lines = []
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

def parse_uploaded_file(file) -> List[Tuple[str,str]]:
    """
    Supports:
      - FASTA/TXT (single or multi-FASTA)
      - GenBank (.gb/.gbk)
      - SnapGene (.dna) if `snapgene-reader` is installed (best-effort; FASTA recommended)
    """
    name = file.name
    lower = name.lower()
    content = file.read()

    if lower.endswith((".fasta",".fa",".fna",".txt")):
        try:
            text = content.decode("utf-8", errors="ignore")
        except Exception:
            text = ""
        recs = parse_fasta_text(text)
        # If it's a single raw sequence (Gene1), use filename-derived name
        if len(recs) == 1 and recs[0][0] == "Gene1":
            recs = [(derive_gene_name(name.rsplit(".",1)[0]), recs[0][1])]
        return recs

    if lower.endswith((".gb",".gbk",".genbank")):
        import io
        bio = io.BytesIO(content)
        records = []
        for rec in SeqIO.parse(bio, "genbank"):
            rec_name = rec.name or rec.id or name.rsplit(".",1)[0]
            records.append((derive_gene_name(rec_name), normalize_seq(str(rec.seq))))
        return records

    if lower.endswith(".dna"):
        try:
            from snapgene_reader import snapgene_file_to_dict  # type: ignore
        except Exception as e:
            raise ValueError(
                "SnapGene .dna 解析可能不稳定。若失败，推荐在 SnapGene 中导出 FASTA 后导入。\n\n"
                "若仍想尝试解析，请安装依赖：\n"
                "  python -m pip install snapgene-reader==0.1.23"
            ) from e
        import io
        bio = io.BytesIO(content)
        d = snapgene_file_to_dict(bio)
        seq = normalize_seq(d.get("seq", "") or "")
        if not seq:
            raise ValueError("未从 .dna 解析到序列。建议导出 FASTA 后导入。")
        rec_name = d.get("name") or name.rsplit(".", 1)[0]
        return [(derive_gene_name(rec_name), seq)]

    raise ValueError(f"Unsupported file type: {name}")
