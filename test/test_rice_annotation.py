import json
import tempfile
import unittest
from pathlib import Path

from rice_annotation import (
    RiceAnnotationResolver,
    RiceAnnotationSchemaError,
    RiceAnnotationSourceError,
)


PRIMARY_HEADER = "geneID\tgene_symbol\tAtBestHit\tArabidopsis_symbol\n"
CURATED_HEADER = "Gene_ID\tGene_symbol\n"


class RiceAnnotationResolverTests(unittest.TestCase):
    def _paths(self, root: Path, primary_rows: str, curated_rows: str = ""):
        primary = root / "primary.tsv"
        curated = root / "curated.tsv"
        primary.write_text(PRIMARY_HEADER + primary_rows, encoding="utf-8")
        curated.write_text(CURATED_HEADER + curated_rows, encoding="utf-8")
        return primary, curated

    def test_official_symbol_has_highest_priority(self):
        with tempfile.TemporaryDirectory() as tmp:
            primary, curated = self._paths(
                Path(tmp),
                "LOC_Os01g00010\tOsPRIMARY; OsAlias\tHY5_AT5G11260\tHY5\n",
                "LOC_Os01g00010\tOsCURATED\n",
            )
            result = RiceAnnotationResolver(primary, curated_path=curated).resolve(
                "LOC_Os01g00010"
            )
        self.assertEqual(result.resolved_symbol, "OsPRIMARY")
        self.assertEqual(result.official_rice_symbol, "OsPRIMARY")
        self.assertEqual(result.fallback_type, "OFFICIAL_RICE_SYMBOL")
        self.assertEqual(result.source, "PRIMARY_OFFICIAL:gene_symbol")
        self.assertFalse(result.warning)

    def test_curated_symbol_precedes_arabidopsis_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            primary, curated = self._paths(
                Path(tmp),
                "LOC_Os01g00010\tunknown\tMYB79_AT4G13480\tMYB79\n",
                "LOC_Os01g00010\tOsCURATED\n",
            )
            result = RiceAnnotationResolver(primary, curated_path=curated).resolve(
                "LOC_Os01g00010"
            )
        self.assertEqual(result.resolved_symbol, "OsCURATED")
        self.assertEqual(result.fallback_type, "CURATED_RICE_SYMBOL")
        self.assertIn("PRIMARY_OFFICIAL_SYMBOL_UNAVAILABLE", result.warning)

    def test_arabidopsis_homolog_uses_shared_os_at_naming(self):
        with tempfile.TemporaryDirectory() as tmp:
            primary, curated = self._paths(
                Path(tmp),
                "LOC_Os01g00010\tunknown\tMYB79_AT4G13480\tMYB79\n",
            )
            result = RiceAnnotationResolver(primary, curated_path=curated).resolve(
                "LOC_Os01g00010"
            )
        self.assertEqual(result.resolved_symbol, "Os_At_MYB79")
        self.assertIsNone(result.official_rice_symbol)
        self.assertEqual(result.at_homolog_id, "AT4G13480")
        self.assertEqual(result.at_homolog_symbol, "MYB79")
        self.assertEqual(result.fallback_type, "ARABIDOPSIS_HOMOLOG")
        self.assertIn("inferred from Arabidopsis homolog", result.warning)

    def test_missing_annotation_retains_stable_gene_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            primary, curated = self._paths(Path(tmp), "")
            result = RiceAnnotationResolver(primary, curated_path=curated).resolve(
                "LOC_Os09g99999"
            )
        self.assertEqual(result.resolved_symbol, "LOC_Os09g99999")
        self.assertEqual(result.fallback_type, "STABLE_GENE_ID")
        self.assertEqual(result.source, "STABLE_GENE_ID")
        self.assertIn("ANNOTATION_UNAVAILABLE", result.warning)

    def test_input_hint_conflict_warns_without_override(self):
        with tempfile.TemporaryDirectory() as tmp:
            primary, curated = self._paths(
                Path(tmp),
                "LOC_Os02g10860\tOsHY5L1; OsbZIP18\tHY5_AT5G11260\tHY5\n",
            )
            result = RiceAnnotationResolver(primary, curated_path=curated).resolve(
                "LOC_Os02g10860", input_hint="OsHY5"
            )
        self.assertEqual(result.resolved_symbol, "OsHY5L1")
        self.assertTrue(result.hint_conflict)
        self.assertIn("INPUT_HINT_CONFLICT", result.warning)

    def test_missing_value_hint_is_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            primary, curated = self._paths(
                Path(tmp),
                "LOC_Os01g00010\tOsPRIMARY\tMYB79_AT4G13480\tMYB79\n",
            )
            result = RiceAnnotationResolver(primary, curated_path=curated).resolve(
                "LOC_Os01g00010", input_hint="unknown"
            )
        self.assertIsNone(result.input_hint)
        self.assertFalse(result.hint_conflict)
        self.assertFalse(result.warning)

    def test_missing_required_column_blocks_loading(self):
        with tempfile.TemporaryDirectory() as tmp:
            primary = Path(tmp) / "primary.tsv"
            primary.write_text("geneID\tgene_symbol\nLOC_Os01g00010\tOsTEST\n", encoding="utf-8")
            with self.assertRaises(RiceAnnotationSchemaError):
                RiceAnnotationResolver(primary, curated_path=None)

    def test_unreadable_primary_source_blocks_loading(self):
        with tempfile.TemporaryDirectory() as tmp:
            missing = Path(tmp) / "missing.tsv"
            with self.assertRaises(RiceAnnotationSourceError):
                RiceAnnotationResolver(missing, curated_path=None)



if __name__ == "__main__":
    unittest.main()
