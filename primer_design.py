import re
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction

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
    if stop in {"TAA","TAG","TGA"}:
        return cds[:-3], True
    return cds, False

def tm_nn(seq: str, Na=50.0, Mg=1.5, dNTPs=0.2, primer_conc=250.0):
    # Tm_NN expects Na in mM
    return mt.Tm_NN(seq, Na=Na, Mg=Mg, dNTPs=dNTPs, dnac1=primer_conc, dnac2=primer_conc)

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
    Na: float = 50.0,
    Mg: float = 1.5,
    dNTPs: float = 0.2,
    primer_conc: float = 250.0,
):
    seq = normalize_seq(seq)
    warnings = []

    if remove_stop:
        seq2, trimmed = trim_stop_codon(seq)
        if trimmed:
            warnings.append("Stop codon removed")
        seq = seq2

    if len(seq) < max_len * 2:
        warnings.append("Sequence is short; primer design may be suboptimal")

    homo_f = normalize_seq(homo_f)
    homo_r = normalize_seq(homo_r)

    def pick_forward():
        best = None
        for L in range(min_len, max_len + 1):
            s = seq[:L]
            if len(s) < L:
                continue
            t = tm_nn(s, Na=Na, Mg=Mg, dNTPs=dNTPs, primer_conc=primer_conc)
            score = abs(t - target_tm)
            if best is None or score < best[0]:
                best = (score, s, t)
        return best

    def pick_reverse():
        best = None
        for L in range(min_len, max_len + 1):
            s = seq[-L:]
            if len(s) < L:
                continue
            s_rc = revcomp(s)
            t = tm_nn(s_rc, Na=Na, Mg=Mg, dNTPs=dNTPs, primer_conc=primer_conc)
            score = abs(t - target_tm)
            if best is None or score < best[0]:
                best = (score, s_rc, t)
        return best

    fb = pick_forward()
    rb = pick_reverse()

    gene_bind_f = fb[1]
    gene_bind_r = rb[1]
    tm_f = float(fb[2])
    tm_r = float(rb[2])

    if not (tm_min <= tm_f <= tm_max):
        warnings.append(f"TmF out of range ({tm_f:.2f})")
    if not (tm_min <= tm_r <= tm_max):
        warnings.append(f"TmR out of range ({tm_r:.2f})")

    gc_f = float(gc_fraction(gene_bind_f) * 100.0) if gene_bind_f else None
    gc_r = float(gc_fraction(gene_bind_r) * 100.0) if gene_bind_r else None

    final_f = homo_f + gene_bind_f
    final_r = homo_r + gene_bind_r

    return {
        "gene_name": gene_name,
        "gene_bind_f": gene_bind_f,
        "gene_bind_r": gene_bind_r,
        "final_f": final_f,
        "final_r": final_r,
        "tm_f": tm_f,
        "tm_r": tm_r,
        "gc_f": gc_f,
        "gc_r": gc_r,
        "warnings": "; ".join(warnings),
    }
