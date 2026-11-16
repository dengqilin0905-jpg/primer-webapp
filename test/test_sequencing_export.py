import unittest

from openpyxl import load_workbook
import io

from project_contract import LEGACY_CONFIRMED_SEQUENCING_TEMPLATE
from sequencing_export import (
    DEFAULT_COLONY_COUNT,
    PROTECTED_SHEET,
    UPLOAD_SHEET,
    SequencingExportError,
    export_sequencing_workbook,
    resolve_sequencing_template,
)


def item(label="PHW_OsTF26", reference="AACCGGTTAACCGGTTAA"):
    return {
        "operational_label": label,
        "full_construct_name": "SYNTHETIC_full_construct",
        "reference_primer_name": f"{label}_F",
        "reference_sequence": reference,
        "plasmid_length": 12000,
        "sequencing_length": 1499,
        "resistance": "氨苄",
        "source_map": "/maps/full_construct.gbk",
    }


class SequencingExportTests(unittest.TestCase):
    def test_registered_template_resolver(self):
        resolution = resolve_sequencing_template()
        self.assertEqual(resolution.source_path, str(LEGACY_CONFIRMED_SEQUENCING_TEMPLATE.resolve()))
        self.assertEqual(resolution.source_kind, "global")
        self.assertEqual(len(resolution.sha256), 64)

    def test_default_three_colonies_and_operational_names(self):
        data, report = export_sequencing_workbook([item(), item("PHW_OsRTF16")])
        self.assertEqual(DEFAULT_COLONY_COUNT, 3)
        self.assertEqual(report["sample_count"], 6)
        self.assertEqual(
            report["sample_names"],
            [
                "PHW_OsTF26-1", "PHW_OsTF26-2", "PHW_OsTF26-3",
                "PHW_OsRTF16-1", "PHW_OsRTF16-2", "PHW_OsRTF16-3",
            ],
        )
        self.assertNotIn(b"SYNTHETIC_full_construct-1", data)

    def test_reference_sequence_is_preserved_and_named_independently(self):
        reference = "AACCGGTTAACCGGTTAA"
        data, report = export_sequencing_workbook([item(reference=reference)])
        workbook = load_workbook(io.BytesIO(data), data_only=False)
        self.assertEqual(workbook[UPLOAD_SHEET].cell(4, 4).value, reference)
        self.assertEqual(report["reference_primers"], ["PHW_OsTF26_F"])

    def test_explicit_colony_override(self):
        _, report = export_sequencing_workbook([item()], colony_count=2)
        self.assertEqual(report["sample_names"], ["PHW_OsTF26-1", "PHW_OsTF26-2"])

    def test_invalid_reference_name_or_sequence_fails(self):
        wrong_name = item()
        wrong_name["reference_primer_name"] = "PHW_TF26_F"
        with self.assertRaisesRegex(SequencingExportError, "reference_primer_name"):
            export_sequencing_workbook([wrong_name])
        wrong_sequence = item(reference="ACGT-N")
        with self.assertRaisesRegex(SequencingExportError, "A/C/G/T"):
            export_sequencing_workbook([wrong_sequence])

    def test_duplicate_operational_label_fails(self):
        with self.assertRaisesRegex(SequencingExportError, "duplicate operational_label"):
            export_sequencing_workbook([item(), item()])

    def test_stale_template_rows_are_cleared(self):
        first,_ = export_sequencing_workbook([item('SYN_OsDemo1'),item('SYN_OsDemo2')])
        data, report = export_sequencing_workbook([item()],template_bytes=first)
        workbook = load_workbook(io.BytesIO(data), data_only=False)
        sheet = workbook[UPLOAD_SHEET]
        observed = [
            [sheet.cell(row, column).value for column in range(1, 11)]
            for row in range(4, sheet.max_row + 1)
            if any(sheet.cell(row, column).value is not None for column in range(1, 11))
        ]
        self.assertEqual(len(observed), 3)
        self.assertEqual(report["sample_count"], 3)

    def test_template_structure_and_validations_are_preserved(self):
        data, report = export_sequencing_workbook([item()])
        workbook = load_workbook(io.BytesIO(data), data_only=False)
        self.assertEqual(workbook[PROTECTED_SHEET].sheet_state, "hidden")
        self.assertEqual(len(workbook.sheetnames), 2)
        self.assertEqual(report["validation_count"], 5)


if __name__ == "__main__":
    unittest.main()
