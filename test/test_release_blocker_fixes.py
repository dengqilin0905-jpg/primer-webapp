"""Focused blocker regressions; synthetic sequences only."""
import unittest,io,json,zipfile,struct
from pathlib import Path
from unittest.mock import patch
from Bio.Seq import Seq
from openpyxl import load_workbook
from test_sequencing_verification import fixture,event
from test_pre_release_acceptance import fusion_fixture,supported_event
from test_sequencing_export import item
from verification_import import bam_metrics
from verification_align import verify_event
from verification_export import project_exports
from verification_html import alignment_columns
from sequencing_export import export_sequencing_workbook,validate_sequencing_workbook,SequencingExportError,UPLOAD_SHEET

class ScientificBlockerTests(unittest.TestCase):
    def deleted_event(self,depth=12,blocks=True,reverse=True):
        ref=fixture();q=ref.sequence[482:600]+ref.sequence[630:758];e=event(q,depth=depth)
        if blocks:e.bam_metrics['read_blocks']=[{'start0':0,'end0':len(q),'sequence':q,'read_name':'read'+str(i),'reverse':bool(i%2) if reverse else False} for i in range(6)]
        return ref,e
    def test_empty_evidence_cannot_confirm_deletion(self):
        ref,e=self.deleted_event(blocks=False);x=verify_event(e,ref,confirmed=True)
        self.assertEqual(x['target_status'],'INCOMPLETE')
        self.assertFalse(next(d for d in x['differences'] if d['type']=='DEL')['support_evidence']['supported'])
    def test_one_x_internal_30bp_deletion_incomplete(self):
        ref,e=self.deleted_event(depth=1);self.assertEqual(verify_event(e,ref,confirmed=True)['target_status'],'INCOMPLETE')
    def test_actual_bidirectional_spanning_reads_confirm_key_deletion(self):
        ref,e=self.deleted_event();x=verify_event(e,ref,confirmed=True)
        self.assertEqual(x['target_status'],'FAIL');self.assertTrue(x['expected_design_mismatch_confirmed'])
        d=next(d for d in x['differences'] if d['type']=='DEL');self.assertEqual(d['support_evidence']['unique_spanning_reads'],6)
        self.assertEqual(x['actual_sample_identity_status'],'NOT_ESTABLISHED')
    def test_single_direction_cannot_confirm_deletion(self):
        ref,e=self.deleted_event(reverse=False);self.assertEqual(verify_event(e,ref,confirmed=True)['target_status'],'INCOMPLETE')
    def test_unconfirmed_reference_default_is_review(self):
        ref,e=self.deleted_event();x=verify_event(e,ref)
        self.assertEqual(x['target_status'],'REVIEW');self.assertFalse(x['reference_assignment_confirmed'])
    def test_wrong_read_anchor_sequence_is_not_support(self):
        ref,e=self.deleted_event()
        for b in e.bam_metrics['read_blocks']:b['sequence']='A'*len(e.consensus)
        self.assertEqual(verify_event(e,ref,confirmed=True)['target_status'],'INCOMPLETE')

    def test_reverse_deletion_uses_original_query_anchors(self):
        ref,e=self.deleted_event();q=str(Seq(e.consensus).reverse_complement());e=event(q)
        e.bam_metrics['read_blocks']=[{'start0':0,'end0':len(q),'sequence':q,'read_name':'read'+str(i),'reverse':bool(i%2)} for i in range(6)]
        x=verify_event(e,ref,confirmed=True)
        self.assertEqual(x['orientation'],'-');self.assertEqual(x['target_status'],'FAIL')
    def test_binary_bam_cigar_blocks_decode_actual_read_bases(self):
        header=b'@SQ\tSN:sample\tLN:40\n';data=b'BAM\x01'+struct.pack('<i',len(header))+header+struct.pack('<i',1)+struct.pack('<i',7)+b'sample\0'+struct.pack('<i',40)
        for i in range(2):
            name=('read'+str(i)).encode()+b'\0';sequence='TT'+'ACGT'*4+'G'+'TGCA'*4
            cigar=b''.join(struct.pack('<I',(n<<4)|op) for n,op in [(2,4),(16,0),(1,1),(16,0)])
            alphabet='=ACMGRSVTWYHKDBN';packed=bytes((alphabet.index(sequence[j])<<4)|(alphabet.index(sequence[j+1]) if j+1<len(sequence) else 0) for j in range(0,len(sequence),2))
            core=struct.pack('<iiIIiiii',0,4,(60<<8)|len(name),((16 if i else 0)<<16)|4,len(sequence),-1,-1,0)
            record=core+name+cigar+packed+bytes([40])*len(sequence);data+=struct.pack('<i',len(record))+record
        metrics=bam_metrics(data);blocks=metrics['read_blocks']
        self.assertEqual([b['sequence'] for b in blocks],['ACGT'*4,'TGCA'*4]*2)
        self.assertEqual([(b['start0'],b['end0']) for b in blocks],[(4,20),(20,36)]*2)
        self.assertEqual([b['reverse'] for b in blocks],[False,False,True,True])

class AntibioticBlockerTests(unittest.TestCase):
    def test_explicit_amp_alias_maps_to_template_value(self):
        data,report=export_sequencing_workbook([{**item(),'resistance':'AmpR'}]);wb=load_workbook(io.BytesIO(data))
        self.assertEqual(wb[UPLOAD_SHEET]['G4'].value,'氨苄');self.assertEqual(report['sample_count'],3)
    def test_unknown_value_blocked(self):
        with self.assertRaises(SequencingExportError):export_sequencing_workbook([{**item(),'resistance':'unknown_DRUG'}])
    def test_illegal_output_cell_fails_readback(self):
        data,_=export_sequencing_workbook([item()]);wb=load_workbook(io.BytesIO(data));wb[UPLOAD_SHEET]['G4']='AmpR';b=io.BytesIO();wb.save(b)
        with self.assertRaises(SequencingExportError):validate_sequencing_workbook(b.getvalue(),[item()])
    def test_validation_covers_all_expanded_rows(self):
        data,_=export_sequencing_workbook([item('PHW_OsTF'+str(i)) for i in range(30)]);wb=load_workbook(io.BytesIO(data));ws=wb[UPLOAD_SHEET]
        self.assertTrue(any(str(v.formula1).lstrip('=')=='Antibiotic' and 'G93' in v for v in ws.data_validations.dataValidation))

    def test_alias_blocked_when_target_not_allowed_by_template(self):
        from sequencing_export import resolve_sequencing_template,antibiotic_options
        wb=load_workbook(io.BytesIO(resolve_sequencing_template().template_bytes))
        for sheet,area in wb.defined_names['Antibiotic'].destinations:
            for row in wb[sheet][area]:
                for cell in row:
                    if cell.value=='氨苄':cell.value='different_official_value'
        self.assertNotIn('氨苄',antibiotic_options(wb));b=io.BytesIO();wb.save(b)
        with self.assertRaises(SequencingExportError):export_sequencing_workbook([{**item(),'resistance':'AmpR'}],template_bytes=b.getvalue())

class HtmlBlockerTests(unittest.TestCase):
    def test_real_bases_gaps_coordinates_and_features(self):
        ref=fixture();raw=ref.sequence[482:600]+'A'+ref.sequence[600:758];e=event(raw);x=verify_event(e,ref,confirmed=True)
        cols=alignment_columns(x,ref,e)
        self.assertEqual(''.join(c['query'] for c in cols if c['query']!='-'),raw)
        self.assertTrue(any(c['kind']=='insertion' and c['ref']=='-' for c in cols))
        html=project_exports([x],[ref],[e])['html'].decode()
        for text in ['base-alignment','data-ref0','data-query0','insertion','synthetic non-CDS insert','Reference +','Query + original']:
            self.assertIn(text,html)
    def test_reverse_coordinates_correspond_to_original_fasta(self):
        ref=fixture();raw=str(Seq(ref.sequence[482:758]).reverse_complement());e=event(raw);x=verify_event(e,ref,confirmed=True)
        for c in alignment_columns(x,ref,e):
            if c['query_pos0'] is not None:self.assertEqual(str(Seq(raw[c['query_pos0']]).complement()),c['query'])
        for seg in x['segments']:
            self.assertEqual(seg['query_start0'],len(raw)-seg['query_oriented_end0'])
    def test_deletion_is_visible_with_actual_context(self):
        ref=fixture();q=ref.sequence[482:600]+ref.sequence[630:758];e=event(q);x=verify_event(e,ref,confirmed=True)
        cols=alignment_columns(x,ref,e);self.assertEqual(sum(c['kind']=='deletion' for c in cols),30)
        html=project_exports([x],[ref],[e])['html'].decode();self.assertIn('Difference contexts',html);self.assertIn('Deletion support',html)

    def test_mismatch_highlights_actual_reference_and_query_bases(self):
        ref=fixture();q=list(ref.sequence[482:758]);q[100]='A' if q[100]!='A' else 'C';e=event(''.join(q));x=verify_event(e,ref,confirmed=True)
        mismatch=[c for c in alignment_columns(x,ref,e) if c['kind']=='mismatch']
        self.assertEqual(len(mismatch),1);self.assertNotEqual(mismatch[0]['ref'],mismatch[0]['query'])
        html=project_exports([x],[ref],[e])['html'].decode();self.assertIn('class="base mismatch"',html)

class HistoryFilterTests(unittest.TestCase):
    def test_batch_construct_colony_filters_existing_history(self):
        from streamlit.testing.v1 import AppTest
        app=Path(__file__).resolve().parents[1]/'app'
        rows=[{'company_id':'S1','sample_name':'X-1','batch':'batch-a','construct':'X','colony_id':'colony-A','colony_suffix':'1','target_status':'REVIEW','whole_plasmid_status':'INCOMPLETE','reference_sha256':'a'*64,'policy_sha256':'b'*64,'regions':[],'differences':[]},
              {'company_id':'S2','sample_name':'Y-2','batch':'batch-b','construct':'Y','colony_id':'colony-B','colony_suffix':'2','target_status':'INCOMPLETE','whole_plasmid_status':'INCOMPLETE','reference_sha256':'a'*64,'policy_sha256':'b'*64,'regions':[],'differences':[]}]
        script="import sys,streamlit as st\nsys.path.insert(0,"+repr(str(app))+")\nfrom verification_ui import sequencing_verification_page\nst.session_state.setdefault('verification_987654_history_results',"+repr(rows)+")\nsequencing_verification_page({'id':987654},None)"
        at=AppTest.from_string(script,default_timeout=20).run();self.assertFalse(at.exception)
        def shown():return next(df.value for df in at.dataframe if 'company_id' in df.value.columns)
        self.assertEqual(len(shown()),2)
        at.selectbox(key='verification_987654_filter_batch').select('batch-a').run();self.assertEqual(shown()['company_id'].tolist(),['S1'])
        at.selectbox(key='verification_987654_filter_batch').select('全部').run()
        at.selectbox(key='verification_987654_filter_construct').select('Y').run();self.assertEqual(shown()['company_id'].tolist(),['S2'])
        at.text_input(key='verification_987654_filter_colony').input('colony-B').run();self.assertEqual(len(shown()),1)
        at.text_input(key='verification_987654_filter_colony').input('1').run();self.assertFalse(any('company_id' in df.value.columns for df in at.dataframe))
        self.assertFalse(at.exception)

class ContractStateTests(unittest.TestCase):
    def test_rerun_persists_and_changed_inputs_invalidate_downloads(self):
        from streamlit.testing.v1 import AppTest
        app=Path(__file__).resolve().parents[1]/'app'
        script="""import sys,io,zipfile
from unittest.mock import patch
import streamlit as st
sys.path.insert(0,APP)
from project_ui import project_contract_page
st.session_state.setdefault('payload',b'{}')
b=io.BytesIO()
with zipfile.ZipFile(b,'w') as z:z.writestr('results/sequencing/sequencing_order.xlsx',b'xlsx-fixture')
r={'summary':{'project':'synthetic'},'conformance':{'status':'PASS'},'normalized_canonical_result_sha256':'synthetic'}
with patch('project_ui.st.file_uploader',return_value=io.BytesIO(st.session_state['payload'])),patch('project_ui.build_project_package',return_value=(b.getvalue(),r)):
    project_contract_page({'id':987654},None)
""".replace('APP',repr(str(app)))
        at=AppTest.from_string(script,default_timeout=20).run();at.button[0].click().run();self.assertFalse(at.exception);self.assertEqual(len(at.get('download_button')),2)
        at.run();self.assertEqual(len(at.get('download_button')),2)
        at.session_state['payload']=b'changed';at.run();self.assertEqual(len(at.get('download_button')),0)

if __name__=='__main__':unittest.main()
