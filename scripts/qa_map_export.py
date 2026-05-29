import io
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "app"
DEPS = APP / "_deps"
for p in [str(DEPS), str(APP)]:
    if p not in sys.path:
        sys.path.insert(0, p)

from Bio import SeqIO

from excel_export import TEMPLATE_PATH, export_from_template_bytes, validate_primer_rows
from genbank_export import build_genbank_files
from primer_design import design_primers_for_gene, normalize_seq
from seq_io import parse_uploaded_records
from snapgene_export import build_snapgene_dna_bundle, find_snapgene_cli


class NamedBytesIO(io.BytesIO):
    def __init__(self, path: Path):
        super().__init__(path.read_bytes())
        self.name = path.name


def _load_one(path: Path) -> dict:
    records = parse_uploaded_records(NamedBytesIO(path))
    if not records:
        raise AssertionError(f"No sequence parsed from {path}")
    rec = records[0]
    rec["name"] = path.stem
    return rec


def _feature_labels(record) -> set[str]:
    labels = set()
    for feature in record.features:
        for key in ["label", "gene", "note"]:
            for value in feature.qualifiers.get(key, []):
                labels.add(str(value))
    return labels


def _assert_feature(record, ftype: str | None, label_text: str):
    for feature in record.features:
        labels = []
        for key in ["label", "gene", "note"]:
            labels.extend(str(x) for x in feature.qualifiers.get(key, []))
        if (ftype is None or feature.type == ftype) and any(label_text in x for x in labels):
            return
    labels = sorted(_feature_labels(record))
    raise AssertionError(f"Missing feature {ftype or '*'} / {label_text}. Seen labels: {labels[:30]}")


def _build_row(backbone: dict, gene: dict) -> dict:
    out = design_primers_for_gene(
        gene_name=gene["name"],
        seq=gene["sequence"],
        homo_f="CGCTCTAGAACTAGTGGATCC",
        homo_r="GATAAGCTTGATATCGAATTC",
        seq_type="CDS",
        remove_stop=False,
        min_len=18,
        max_len=30,
        tm_min=55.0,
        tm_max=65.0,
        target_tm=60.0,
        gc_min=40.0,
        gc_max=60.0,
        restriction_site_f="GGATCC",
        restriction_site_r="GAATTC",
        junction_extra_bases_f=0,
        junction_extra_bases_r=0,
        cds_recovery_enable=True,
    )
    out.update(
        {
            "primer_name_f": f"pGreenII-62SK_{gene['name']}-F",
            "primer_name_r": f"pGreenII-62SK_{gene['name']}-R",
            "backbone_name": "pGreenII-62SK",
            "backbone_label": "pGreenII-62SK [BamHI/EcoRI]",
            "gene_name": gene["name"],
            "seq_type": "CDS",
            "remove_stop": False,
            "construct_type": "CDS_no_tag",
            "homo_f": "CGCTCTAGAACTAGTGGATCC",
            "homo_r": "GATAAGCTTGATATCGAATTC",
            "vector_sequence": backbone["sequence"],
            "insertion_mode": "restriction_site_join",
            "enzyme_sites": "BamHI/EcoRI",
            "restriction_site_f": "GGATCC",
            "restriction_site_r": "GAATTC",
            "junction_extra_bases_f": 0,
            "junction_extra_bases_r": 0,
            "vector_map_note": "QA map export",
            "vector_features_json": json.dumps(backbone.get("features") or [], ensure_ascii=False),
            "vector_topology": backbone.get("topology") or "circular",
            "vector_source_format": backbone.get("source_format") or "snapgene",
            "vector_record_name": backbone.get("record_name") or "",
            "insert_feature_type": "CDS",
            "insert_feature_name": gene["name"],
            "homology_feature_type": "misc_feature",
        }
    )
    return out


def main():
    out_dir = ROOT / ".codex-qa" / "map_export"
    out_dir.mkdir(parents=True, exist_ok=True)

    backbone = _load_one(ROOT / "test" / "pGreenII-62SK.dna")
    genes = [_load_one(ROOT / "test" / "ZmMYBR26.dna"), _load_one(ROOT / "test" / "ZmMYBS1.dna")]
    rows = [_build_row(backbone, gene) for gene in genes]

    validate_primer_rows(rows)
    xlsx = export_from_template_bytes(rows, template_bytes=TEMPLATE_PATH.read_bytes(), profile=None)
    xlsx_path = out_dir / "Sangon_Primers_Template.xlsx"
    xlsx_path.write_bytes(xlsx)
    if xlsx_path.stat().st_size <= 0:
        raise AssertionError("Excel export is empty")

    gbk_files = build_genbank_files(rows)
    for name, data in gbk_files:
        path = out_dir / name
        path.write_bytes(data)
        record = SeqIO.read(io.StringIO(data.decode("utf-8")), "genbank")
        _assert_feature(record, "CDS", Path(name).stem.split("_")[1])
        _assert_feature(record, "misc_feature", "F_homology_arm")
        _assert_feature(record, "misc_feature", "R_homology_arm")
        _assert_feature(record, "misc_feature", "BamHI")
        _assert_feature(record, "misc_feature", "EcoRI")
        _assert_feature(record, "promoter", "CaMV 35S promoter")
        _assert_feature(record, "CDS", "KanR")

    snapgene_exe, detail = find_snapgene_cli()
    if snapgene_exe is None:
        print(f"SNAPGENE_SKIPPED {detail}")
    else:
        try:
            dna_name, dna_data, _, used_exe = build_snapgene_dna_bundle(rows)
            dna_path = out_dir / dna_name
            dna_path.write_bytes(dna_data)
            if dna_path.stat().st_size <= 0:
                raise AssertionError("SnapGene DNA export is empty")
            print(f"SNAPGENE_OK {used_exe} -> {dna_path}")
        except Exception as exc:
            print(f"SNAPGENE_FAILED {snapgene_exe}: {exc}")

    print(f"OK rows={len(rows)} gbk={len(gbk_files)} excel={xlsx_path}")


if __name__ == "__main__":
    main()
