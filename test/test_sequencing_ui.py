import unittest,io,json,random,inspect
from types import SimpleNamespace
from unittest.mock import patch
from Bio.Seq import Seq
from openpyxl import load_workbook
from sequencing_ui import items_from_design
from sequencing_export import export_sequencing_workbook
from identity_audit import IdentityEvidence

def row():
    rng=random.Random(31);v=''.join(rng.choice('ACGT') for _ in range(1800));a,b=500,560
    return {'gene_name':'OsSynthetic','backbone_name':'PHW_3xHA','vector_sequence':v,
            'homo_f':v[a-18:a],'homo_r':str(Seq(v[b:b+18]).reverse_complement()),
            'insert_seq':'ATG'+'GCT'*59+'TAA','vector_features_json':'[]','vector_topology':'circular'}

def meta():return {'symbol':'OsSynthetic','species':'Os','canonical_identity':'','resistance':'AmpR'}

class SequencingAdaptersTests(unittest.TestCase):
    def test_existing_design_export_preserves_sequence_and_three_colonies(self):
        original=row();before=json.dumps(original,sort_keys=True)
        items=items_from_design([original],[meta()]);data,report=export_sequencing_workbook(items)
        self.assertEqual(report['sample_count'],3);self.assertEqual(json.dumps(original,sort_keys=True),before)
        self.assertEqual(items[0]['reference_sequence'],original['homo_f'])
        self.assertEqual(items[0]['sequencing_length'],len(original['insert_seq'])+18)
        wb=load_workbook(io.BytesIO(data));ws=wb.worksheets[0]
        self.assertEqual([ws.cell(i,1).value for i in range(4,7)],['PHW_OsSynthetic-1','PHW_OsSynthetic-2','PHW_OsSynthetic-3'])
    def test_legacy_21bp_arm_not_silently_truncated(self):
        data=row();data['homo_f']='A'*21
        with self.assertRaises(ValueError):items_from_design([data],[meta()])
    def test_numbered_identity_requires_matching_registry_locus(self):
        m={**meta(),'canonical_identity':'TF999','gene_id':'LOC_SYNTHETIC'}
        evidence=[IdentityEvidence('LOC_OTHER','TF999','identity_registry','synthetic')]
        with patch('sequencing_ui.load_identity_tsv',return_value=evidence):
            with self.assertRaises(ValueError):items_from_design([row()],[m])
        evidence=[IdentityEvidence('LOC_SYNTHETIC','TF999','identity_registry','synthetic')]
        with patch('sequencing_ui.load_identity_tsv',return_value=evidence):
            self.assertEqual(items_from_design([row()],[m])[0]['operational_label'],'PHW_OsTF999')
    def test_resistance_and_species_must_be_explicit(self):
        with self.assertRaises(ValueError):items_from_design([row()],[{**meta(),'resistance':''}])
    def test_page_connections_preserve_contract_zip(self):
        from pathlib import Path
        app=Path(__file__).resolve().parents[1]/'app'
        self.assertIn('sequencing_download(all_result_rows',(app/'design_ui.py').read_text())
        self.assertIn('sequencing_download(primer_rows',(app/'history_ui.py').read_text())
        s=(app/'project_ui.py').read_text();self.assertIn('Download canonical project ZIP',s);self.assertIn('Download sequencing Excel',s)
        self.assertIn('Sequencing Verification',(app/'app.py').read_text())

if __name__=='__main__':unittest.main()
