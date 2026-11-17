import json
import tempfile
import unittest
from pathlib import Path

from project_contract import (
    CANONICAL_DIRECTORIES,
    CANONICAL_OUTPUTS,
    canonical_project_paths,
    canonical_species_abbreviation,
    minimal_backbone_name,
)
from project_manifest import ProjectManifestError, load_project_manifest, validate_project_manifest
from project_validator import audit_project


def valid_manifest(project_id="contract_test"):
    return {
        "contract_version": 1,
        "project_id": project_id,
        "project_type": "primer_design",
        "created_at": "2026-09-16T10:00:00+08:00",
        "species": ["Rice"],
        "constructs": [
            {
                "gene_id": "Os01g0000000",
                "symbol": "OsTEST1",
                "canonical_identity": "TF999",
                "operational_label": "PHW_OsTF999",
                "full_construct_name": "PHW_TF999_OsTEST1_3xHA",
                "backbone": "PHW",
                "preset_id": "PHW_3xHA",
                "preset_version": 1,
                "remove_stop": True,
                "physical_plasmid_status": "not_found",
                "overrides": [],
                "status": "planned",
            }
        ],
        "outputs": dict(CANONICAL_OUTPUTS),
        "provenance": {
            "identity_sources": ["identity_master.tsv"],
            "vector_hash": "sha256:test",
            "template_source": "/canonical/order_primer.xlsx",
        },
    }


class ProjectContractTests(unittest.TestCase):
    def test_canonical_paths_are_stable(self):
        paths = canonical_project_paths("/tmp/example")
        self.assertEqual(str(paths["primer_master"]), "/tmp/example/results/primers/primer_master.tsv")
        self.assertEqual(str(paths["order"]), "/tmp/example/results/order/primer_order.xlsx")
        self.assertEqual(str(paths["sequencing"]), "/tmp/example/results/sequencing/sequencing_order.xlsx")

    def test_species_and_backbone_registry(self):
        self.assertEqual(canonical_species_abbreviation("Rice"), "Os")
        self.assertEqual(canonical_species_abbreviation("Nicotiana benthamiana"), "Nb")
        self.assertEqual(minimal_backbone_name("PHW_3xHA"), "PHW")
        self.assertEqual(minimal_backbone_name("pGreenII-62SK"), "62SK")
        with self.assertRaises(ValueError):
            canonical_species_abbreviation("unknown species")

    def test_manifest_loader_and_override_contract(self):
        manifest = valid_manifest()
        self.assertEqual(validate_project_manifest(manifest)["project_id"], "contract_test")
        manifest["constructs"][0]["overrides"] = [
            {
                "override_type": "scientific_qc",
                "override_reason": "reporter boundary requires an additional check",
                "override_value": "check_promoter_boundary",
            }
        ]
        self.assertEqual(len(validate_project_manifest(manifest)["constructs"][0]["overrides"]), 1)
        manifest["constructs"][0]["overrides"][0]["override_type"] = "silent_behavior_change"
        with self.assertRaises(ProjectManifestError):
            validate_project_manifest(manifest)

    def test_manifest_requires_canonical_outputs(self):
        manifest = valid_manifest()
        manifest["outputs"]["order"] = "results/order/final_v2.xlsx"
        with self.assertRaisesRegex(ProjectManifestError, "outputs.order"):
            validate_project_manifest(manifest)

    def test_manifest_requires_species_aware_operational_label(self):
        manifest = valid_manifest()
        manifest["constructs"][0]["operational_label"] = "PHW_TF999"
        with self.assertRaisesRegex(ProjectManifestError, "PHW_OsTF999"):
            validate_project_manifest(manifest)

    def test_read_only_project_audit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "contract_test"
            root.mkdir()
            for relative in CANONICAL_DIRECTORIES:
                (root / relative).mkdir(parents=True, exist_ok=True)
            for relative in (
                "WORK_LOG.md",
                "construct_plan.tsv",
                "results/primers/primer_master.tsv",
                "results/qc/primer_qc.tsv",
                "results/order/primer_order.xlsx",
                "results/sequencing/sequencing_order.xlsx",
            ):
                (root / relative).write_bytes(b"test")
            (root / "project.yaml").write_text(json.dumps(valid_manifest()), encoding="utf-8")
            before = sorted(str(path.relative_to(root)) for path in root.rglob("*"))
            report = audit_project(root)
            after = sorted(str(path.relative_to(root)) for path in root.rglob("*"))
            self.assertEqual(report["status"], "FAIL")
            self.assertTrue(
                any(
                    item["check"] == "file:inputs/runtime_inputs.json"
                    and item["status"] == "FAIL"
                    for item in report["checks"]
                )
            )
            self.assertEqual(before, after)

            (root / "results/order/primer_order.xlsx").unlink()
            report = audit_project(root)
            self.assertEqual(report["status"], "FAIL")
            self.assertTrue(any(item["check"] == "file:results/order/primer_order.xlsx" for item in report["checks"]))

    def test_loader_reads_json_compatible_yaml(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "project.yaml"
            path.write_text(json.dumps(valid_manifest()), encoding="utf-8")
            self.assertEqual(load_project_manifest(path)["contract_version"], 1)


if __name__ == "__main__":
    unittest.main()
