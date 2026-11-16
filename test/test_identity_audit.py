import tempfile
import unittest
from pathlib import Path

from identity_audit import (
    IdentityAuditError,
    IdentityEvidence,
    audit_identity,
    discover_active_project_evidence,
    identity_key,
    load_identity_tsv,
    next_available_tf_id,
    reserved_identities,
    resolve_identity_column,
)
from project_contract import IDENTITY_PREFLIGHT_SOURCE_KINDS


def evidence_for(
    gene_id,
    identity,
    source_kind,
    physical="unknown",
):
    return IdentityEvidence(gene_id, identity, source_kind, f"/{source_kind}.tsv", physical)


def complete_sources(gene_id="LOC_Os01g00010", identity="OsTF26", physical="uncertain"):
    return [
        evidence_for(gene_id, identity, source_kind, physical if source_kind == "physical_inventory" else "unknown")
        for source_kind in IDENTITY_PREFLIGHT_SOURCE_KINDS
    ]


class IdentityAuditTests(unittest.TestCase):

    def test_prefixed_and_unprefixed_identity_are_equivalent(self):
        self.assertEqual(identity_key("TF26"), "TF26")
        self.assertEqual(identity_key("OsTF26"), "TF26")
        self.assertEqual(identity_key("RTF16"), "RTF16")

    def test_existing_available_plasmid_skips_only_when_requested(self):
        result = audit_identity(
            "LOC_Os01g00010",
            complete_sources(physical="AVAILABLE"),
            requested_identity="TF26", respect_existing_physical=True,
        )
        self.assertEqual(result["canonical_identity"], "OsTF26")
        self.assertEqual(result["decision"], "SKIP_REBUILD")

    def test_missing_plasmid_rebuilds_same_identity(self):
        result = audit_identity(
            "LOC_Os01g00010",
            complete_sources(identity="TF50", physical="FAILED"),
        )
        self.assertEqual(result["canonical_identity"], "TF50")
        self.assertEqual(result["decision"], "REBUILD_EXISTING_IDENTITY")

    def test_uncertain_physical_state_is_nonblocking_warning(self):
        result = audit_identity("LOC_Os01g00010", complete_sources(physical="UNCERTAIN"))
        self.assertEqual(result["decision"], "REBUILD_EXISTING_IDENTITY")
        self.assertEqual(result["physical_plasmid_status"], "unknown")
        self.assertEqual(result["source_status"]["physical_inventory"], "AVAILABLE")

    def test_design_records_do_not_prove_physical_plasmid(self):
        rows = [
            evidence_for("LOC_Os01g00010", "TF50", source_kind)
            for source_kind in IDENTITY_PREFLIGHT_SOURCE_KINDS
        ]
        result = audit_identity("LOC_Os01g00010", rows)
        self.assertEqual(result["decision"], "REBUILD_EXISTING_IDENTITY")

    def test_conflicting_requested_identity_is_blocked(self):
        result = audit_identity(
            "LOC_Os01g00010",
            complete_sources(identity="OsTF26", physical="missing"),
            requested_identity="TF66",
        )
        self.assertEqual(result["decision"], "BLOCK_IDENTITY_CONFLICT")

    def test_missing_hard_master_blocks(self):
        rows = [evidence_for("LOC_OTHER", "TF69", "active_construct_plan")]
        result = audit_identity("LOC_NEW", rows)
        self.assertEqual(result["decision"], "BLOCK_INCOMPLETE_SOURCES")
        self.assertIn("identity_registry", result["missing_source_kinds"])

    def test_unreadable_hard_master_raises(self):
        with self.assertRaises(IdentityAuditError):
            load_identity_tsv("/definitely/missing/missing_identity.tsv", source_kind="identity_registry")

    def test_missing_supplementary_sources_do_not_block(self):
        rows = [evidence_for("LOC_OTHER", "TF69", "identity_registry")]
        result = audit_identity("LOC_NEW", rows)
        self.assertEqual(result["decision"], "READY_FOR_NEW_IDENTITY")
        self.assertEqual(result["source_status"]["active_construct_plan"], "NOT_FOUND")
        self.assertEqual(result["source_status"]["historical_mapping"], "NOT_FOUND")
        self.assertEqual(result["source_status"]["physical_inventory"], "NOT_VERIFIED")

    def test_next_tf_number_uses_registry_and_active_plans_but_not_rtf(self):
        rows = [
            evidence_for("G1", "TF69", "identity_registry"),
            evidence_for("G2", "OsTF27", "active_construct_plan"),
            evidence_for("G3", "RTF279", "historical_mapping"),
        ]
        self.assertEqual(next_available_tf_id(rows), "TF70")

    def _load_one_alias(self, alias, identity="OsTF26"):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        path = Path(tmp.name) / "construct_plan.tsv"
        path.write_text(
            f"gene_id\t{alias}\tphysical_plasmid_status\n"
            f"LOC_Os01g00010\t{identity}\tUNCERTAIN\n",
            encoding="utf-8",
        )
        return load_identity_tsv(path, source_kind="active_construct_plan")

    def test_a_current_uppercase_id_is_readable(self):
        rows = self._load_one_alias("current_ID")
        self.assertEqual(rows[0].canonical_identity, "OsTF26")

    def test_b_current_lowercase_id_is_readable(self):
        rows = self._load_one_alias("current_id")
        self.assertEqual(rows[0].canonical_identity, "OsTF26")

    def test_c_canonical_identity_is_readable(self):
        rows = self._load_one_alias("canonical_identity")
        self.assertEqual(rows[0].canonical_identity, "OsTF26")

    def test_d_multiple_consistent_aliases_use_canonical_priority(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "construct_plan.tsv"
            path.write_text(
                "gene_id\tcurrent_id\tcanonical_identity\n"
                "LOC_Os01g00010\tTF26\tOsTF26\n",
                encoding="utf-8",
            )
            rows = load_identity_tsv(path, source_kind="active_construct_plan")
            self.assertEqual(rows[0].canonical_identity, "OsTF26")

    def test_e_multiple_conflicting_aliases_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "construct_plan.tsv"
            path.write_text(
                "gene_id\tcurrent_id\tcanonical_identity\n"
                "LOC_Os01g00010\tTF26\tTF66\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(IdentityAuditError, "Conflicting identity aliases"):
                load_identity_tsv(path, source_kind="active_construct_plan")

    def test_f_loader_fails_when_schema_is_not_identity_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.tsv"
            path.write_text("sample\tvalue\nA\t1\n", encoding="utf-8")
            with self.assertRaises(IdentityAuditError):
                load_identity_tsv(path, source_kind="identity_registry")

    def test_resolver_is_exact_not_fuzzy(self):
        with self.assertRaises(IdentityAuditError):
            resolve_identity_column(["gene_id", "current_identity_note"])
        with self.assertRaises(IdentityAuditError):
            resolve_identity_column(["gene_ID", "original_ID"])



    def test_active_manifest_is_discovered_without_runtime_plan(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "project"
            active = root / "active_project"
            active.mkdir(parents=True)
            (active / "project.yaml").write_text(
                '{"constructs":[{"gene_id":"LOC_ACTIVE","canonical_identity":"TF70"}]}',
                encoding="utf-8",
            )
            rows = discover_active_project_evidence(projects_root=root)
            self.assertEqual([(row.gene_id, row.canonical_identity) for row in rows], [("LOC_ACTIVE", "TF70")])


if __name__ == "__main__":
    unittest.main()
