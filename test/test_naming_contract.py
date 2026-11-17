import unittest

from naming import (
    build_construct_name,
    build_operational_label,
    build_primer_label,
    build_sequencing_sample_name,
)


class NamingContractTests(unittest.TestCase):
    def test_numbered_identity_adds_species_once(self):
        self.assertEqual(
            build_operational_label("OsMYBS1", species="Rice", identifier="TF26", backbone="PHW_3xHA"),
            "PHW_OsTF26",
        )
        self.assertEqual(
            build_operational_label("OsMYBS1", species="Os", identifier="OsTF26", backbone="PHW"),
            "PHW_OsTF26",
        )
        self.assertEqual(
            build_operational_label("OsMYBS1", species="Rice", identifier="RTF16", backbone="PHW"),
            "PHW_OsRTF16",
        )

    def test_symbol_and_at_homolog_names(self):
        self.assertEqual(
            build_operational_label("MYBS3", species="Rice", backbone="PHW"),
            "PHW_OsMYBS3",
        )
        self.assertEqual(
            build_operational_label("NFXL1", species="Rice", backbone="PHW", at_homolog="NFXL1", official_symbol_available=False),
            "PHW_Os_At_NFXL1",
        )
        self.assertEqual(
            build_operational_label("KUA1", species="Nicotiana benthamiana", backbone="pGreenII-62SK"),
            "62SK_NbKUA1",
        )

    def test_primer_and_sequencing_share_operational_label(self):
        label = build_operational_label("OsMYBS1", species="Rice", identifier="TF26", backbone="PHW")
        self.assertEqual(build_primer_label("OsMYBS1", "F", species="Rice", identifier="TF26", backbone="PHW"), "PHW_OsTF26_F")
        self.assertEqual(build_primer_label("OsMYBS1", "R", species="Rice", identifier="TF26", backbone="PHW"), "PHW_OsTF26_R")
        self.assertEqual([build_sequencing_sample_name(label, i) for i in range(1, 4)], ["PHW_OsTF26-1", "PHW_OsTF26-2", "PHW_OsTF26-3"])

    def test_full_construct_namespace_stays_independent(self):
        full_name = build_construct_name("PHW", "OsMYBS1", identifier="TF26", tag="3xHA")
        operational = build_operational_label("OsMYBS1", species="Rice", identifier="TF26", backbone="PHW")
        self.assertEqual(full_name, "PHW_TF26_OsMYBS1_3xHA")
        self.assertEqual(operational, "PHW_OsTF26")
        self.assertNotEqual(full_name, operational)

    def test_invalid_or_ambiguous_inputs_fail_closed(self):
        with self.assertRaisesRegex(ValueError, "requires species"):
            build_operational_label("MYBS3", identifier="TF26", backbone="PHW")
        with self.assertRaisesRegex(ValueError, "conflicts"):
            build_operational_label("OsMYBS3", species="Arabidopsis", identifier="TF26", backbone="PHW")
        with self.assertRaisesRegex(ValueError, "F or R"):
            build_primer_label("OsMYBS3", "forward", species="Rice", backbone="PHW")
        with self.assertRaisesRegex(ValueError, "positive integer"):
            build_sequencing_sample_name("PHW_OsTF26", 0)


if __name__ == "__main__":
    unittest.main()
