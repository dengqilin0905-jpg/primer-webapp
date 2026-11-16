import unittest,io,json,random,hashlib,zipfile,gzip,struct,tempfile
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
from verification_reference import freeze_reference
from verification_import import Event,parse_zip,merge_events,bam_metrics
from verification_align import verify_event
from verification_export import project_exports
from verification_store import save_verifications,load_events


def fixture(insert=None):
    rng=random.Random(101);original=''.join(rng.choice('ACGT') for _ in range(1800))
    ins=insert or ''.join(rng.choice('ACGT') for _ in range(240));a,b=500,560
    seq=original[:a]+ins+original[b:];rec=SeqRecord(Seq(seq),id='TEST_OsTF999',description='synthetic public fixture')
    rec.annotations={'molecule_type':'DNA','topology':'circular'}
    rec.features=[SeqFeature(FeatureLocation(a,a+len(ins)),type='misc_feature',qualifiers={'label':['synthetic non-CDS insert']})]
    return freeze_reference(rec,backbone_sequence=original,replace_start=a,replace_end=b,insert_sequence=ins,
                            source='synthetic fixture',operational_label='TEST_OsTF999')

def event(query,name='TEST_OsTF999-1',sid='S123',depth=12):
    h=hashlib.sha256((sid+query).encode()).hexdigest()
    return Event(h,sid,name,'synthetic-batch','1',h,h,consensus=query,consensus_sha256=hashlib.sha256(query.encode()).hexdigest(),
                 bam_metrics={'depth':[depth]*len(query),'reference_length':len(query),'forward':6,'reverse':6})

class VerificationTests(unittest.TestCase):
    def test_non_cds_insert_and_junctions_required(self):
        ref=fixture();q=ref.sequence[482:758];x=verify_event(event(q),ref,confirmed=True)
        self.assertEqual(x['target_status'],'PASS');self.assertEqual(x['whole_plasmid_status'],'INCOMPLETE')
        self.assertTrue(any(r['kind']=='insert' for r in x['regions']))
    def test_reverse_and_circular_origin(self):
        ref=fixture();q=ref.sequence[482:758]
        x=verify_event(event(str(Seq(q).reverse_complement())),ref,confirmed=True)
        self.assertEqual(x['target_status'],'PASS');self.assertEqual(x['orientation'],'-')
        q=ref.sequence[-50:]+ref.sequence[:780]
        x=verify_event(event(q),ref,confirmed=True);self.assertEqual(x['target_status'],'PASS')
    def test_missing_junction_and_long_insert_not_pass(self):
        ref=fixture();x=verify_event(event(ref.sequence[500:650]),ref,confirmed=True)
        self.assertEqual(x['target_status'],'INCOMPLETE')
    def test_mutation_in_insert_fails(self):
        ref=fixture();q=list(ref.sequence[482:758]);q[100]='A' if q[100]!='A' else 'C'
        x=verify_event(event(''.join(q)),ref,confirmed=True);self.assertEqual(x['target_status'],'FAIL')
    def test_mixture_never_pass(self):
        ref=fixture();e=event(ref.sequence[482:758]);e.bam_metrics['large_deletions']=[{'start0':80,'length':60,'unique_reads':12,'forward_reads':6,'reverse_reads':6,'fraction':.1}]
        self.assertEqual(verify_event(e,ref,confirmed=True)['target_status'],'REVIEW')
    def test_unconfirmed_and_no_depth(self):
        ref=fixture();e=event(ref.sequence[482:758]);self.assertEqual(verify_event(e,ref)['target_status'],'REVIEW')
        e.bam_metrics={};self.assertEqual(verify_event(e,ref,confirmed=True)['target_status'],'INCOMPLETE')
    def test_images_only(self):
        out=io.BytesIO()
        with zipfile.ZipFile(out,'w') as z:z.writestr('S123/S123-TEST_OsTF999-1.base.png',b'not-rendered-fixture')
        e=parse_zip(out.getvalue())[0];self.assertEqual(verify_event(e,fixture(),confirmed=True)['target_status'],'INCOMPLETE')
    def test_hash_dedup_and_name_collision(self):
        a=event('ACGT',sid='S1');b=event('TGCA',sid='S2');x=merge_events([a],[a,b])
        self.assertEqual(len(x),2);self.assertIn('SAMPLE_NAME_COLLISION',x[0].warnings)
    def test_unsafe_zip(self):
        out=io.BytesIO()
        with zipfile.ZipFile(out,'w') as z:z.writestr('../S123-TEST-1.base.png',b'x')
        with self.assertRaises(ValueError):parse_zip(out.getvalue())
    def test_bad_assembly_rejected(self):
        ref=fixture();rec=__import__('Bio').SeqIO.read(io.StringIO(ref.genbank),'genbank')
        with self.assertRaises(ValueError):freeze_reference(rec,backbone_sequence=ref.sequence,replace_start=0,replace_end=1,insert_sequence='AAAA',source='x',operational_label='TEST_OsTF999')
    def test_private_store_and_feature_exports(self):
        ref=fixture();e=event(ref.sequence[482:758]);x=verify_event(e,ref,confirmed=True)
        outputs=project_exports([x],[ref],[e]);self.assertNotIn(b'<img',outputs['html']);self.assertNotIn(b'<script',outputs['html'])
        with zipfile.ZipFile(io.BytesIO(outputs['zip'])) as z:
            gb=next(n for n in z.namelist() if n.endswith('reference.gbk'))
            self.assertEqual(z.read(gb).decode(),ref.genbank)
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'verification.sqlite3';save_verifications(p,[e],[ref],[x]);save_verifications(p,[e],[ref],[x])
            self.assertEqual(len(load_events(p)),1)
    def test_whole_plasmid_is_independent(self):
        ref=fixture();e=event(ref.sequence);x=verify_event(e,ref,confirmed=True)
        self.assertEqual(x['target_status'],'PASS');self.assertEqual(x['whole_plasmid_status'],'PASS')
    def test_low_depth_cannot_pass(self):
        ref=fixture();self.assertEqual(verify_event(event(ref.sequence[482:758],depth=1),ref,confirmed=True)['target_status'],'INCOMPLETE')


class AdditionalVerificationTests(unittest.TestCase):
    def test_ambiguous_base_is_incomplete_not_mutation(self):
        ref=fixture();q=list(ref.sequence[482:758]);q[100]='N'
        self.assertEqual(verify_event(event(''.join(q)),ref,confirmed=True)['target_status'],'INCOMPLETE')
    def test_repackaged_logically_identical_events_deduplicate(self):
        def zipped(folder):
            b=io.BytesIO()
            with zipfile.ZipFile(b,'w') as z:z.writestr(folder+'/S123-TEST_OsTF999-1.base.png',b'same-content')
            return b.getvalue()
        a=parse_zip(zipped('one'));b=parse_zip(zipped('two'))
        self.assertEqual(len(merge_events(a,b)),1)
    def test_additional_required_element_prevents_pass(self):
        from verification_reference import require_features
        ref=fixture();rec=__import__('Bio').SeqIO.read(io.StringIO(ref.genbank),'genbank')
        rec.features.append(SeqFeature(FeatureLocation(1200,1280),type='promoter',qualifiers={'label':['mandatory promoter']}))
        from dataclasses import replace
        out=io.StringIO();__import__('Bio').SeqIO.write(rec,out,'genbank');ref=replace(ref,genbank=out.getvalue())
        extended=require_features(ref,['mandatory promoter']);e=event(ref.sequence[482:758])
        a=verify_event(e,ref,confirmed=True);b=verify_event(e,extended,confirmed=True)
        self.assertEqual(a['target_status'],'PASS');self.assertEqual(b['target_status'],'INCOMPLETE')
        self.assertNotEqual(a['policy_sha256'],b['policy_sha256'])
        outputs=project_exports([b],[extended],[e])
        with zipfile.ZipFile(io.BytesIO(outputs['zip'])) as z:
            name=next(n for n in z.namelist() if n.endswith('reference_snapshot.json'))
            self.assertIn('mandatory promoter',z.read(name).decode())
    def test_retained_backbone_difference_is_review(self):
        ref=fixture();q=list(ref.sequence[450:810]);q[5]='A' if q[5]!='A' else 'C'
        x=verify_event(event(''.join(q)),ref,confirmed=True)
        self.assertEqual(x['target_status'],'REVIEW');self.assertEqual(x['differences'][0]['classification'],'retained_backbone')
    def test_no_sequence_even_confirmed_reference_is_incomplete(self):
        self.assertEqual(verify_event(event(''),fixture(),confirmed=True)['target_status'],'INCOMPLETE')
    def test_bam_deletion_count_deduplicates_names(self):
        header=b'@SQ\tSN:sample\tLN:120\n';data=b'BAM\x01'+struct.pack('<i',len(header))+header+struct.pack('<i',1)+struct.pack('<i',7)+b'sample\0'+struct.pack('<i',120)
        # Real binary bytes, not the printable escaped header used above.
        header=b'@SQ\tSN:sample\tLN:120\n'.replace(b'\\t',b'\t').replace(b'\\n',b'\n')
        data=b'BAM\x01'.replace(b'\\x01',b'\x01')+struct.pack('<i',len(header))+header+struct.pack('<i',1)+struct.pack('<i',7)+b'sample\0'.replace(b'\\0',b'\0')+struct.pack('<i',120)
        for i in range(8):
            name=('read'+str(i//2)).encode()+b'\0'.replace(b'\\0',b'\0');flag=16 if i%2 else 0
            cigar=b''.join(struct.pack('<I',(n<<4)|op) for n,op in [(30,0),(30,2),(60,0)])
            core=struct.pack('<iiIIiiii',0,0,(60<<8)|len(name),(flag<<16)|3,90,-1,-1,0)
            record=core+name+cigar+bytes([0x11])*45+bytes([40])*90
            data+=struct.pack('<i',len(record))+record
        m=bam_metrics(data);d=m['large_deletions'][0]
        self.assertEqual(d['length'],30);self.assertEqual(d['unique_reads'],4)
    def test_generated_ab1_is_not_independent_quality(self):
        from unittest.mock import patch
        rec=SeqRecord(Seq('ACGT'),id='x');rec.annotations={'abif_raw':{'CMNT1':b'Generated based CWBIO AB1 Creator','RUND1':'2026-01-01'}};rec.letter_annotations['phred_quality']=[50]*4
        b=io.BytesIO()
        with zipfile.ZipFile(b,'w') as z:z.writestr('S123/S123-TEST-1_1.ab1',b'ABI-fixture')
        with patch('verification_import.SeqIO.read',return_value=rec):e=parse_zip(b.getvalue())[0]
        self.assertTrue(e.ab1['generated']);self.assertFalse(e.ab1['independent_quality']);self.assertTrue(e.ab1['constant_q50'])

class FinalPolicyTests(unittest.TestCase):
    def test_per_base_single_strand_gap_cannot_pass(self):
        ref=fixture();e=event(ref.sequence[482:758]);e.bam_metrics['depth_forward']=[6]*len(e.consensus);e.bam_metrics['depth_reverse']=[6]*len(e.consensus)
        e.bam_metrics['depth_reverse'][50]=0
        self.assertEqual(verify_event(e,ref,confirmed=True)['target_status'],'INCOMPLETE')
    def test_short_deletion_subpopulation_cannot_pass(self):
        ref=fixture();e=event(ref.sequence[482:758]);e.bam_metrics['large_deletions']=[{'start0':60,'length':3,'unique_reads':5,'forward_reads':3,'reverse_reads':2,'fraction':.1}]
        self.assertEqual(verify_event(e,ref,confirmed=True)['target_status'],'REVIEW')
    def test_prior_policy_runs_and_snapshots_restore(self):
        from verification_store import load_ledger
        ref=fixture();e=event(ref.sequence[482:758]);a=verify_event(e,ref,confirmed=True)
        b=verify_event(e,ref,confirmed=True,minimum_depth=15)
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'verification.sqlite3';save_verifications(p,[e],[ref],[a]);save_verifications(p,[e],[ref],[b])
            events,refs,results=load_ledger(p)
            self.assertEqual(len(events),1);self.assertEqual(len(results),2);self.assertEqual(refs[0].snapshot(),ref.snapshot())

class WholePlasmidSupportTests(unittest.TestCase):
    def test_full_length_but_low_support_cannot_pass(self):
        ref=fixture();e=event(ref.sequence,depth=1);x=verify_event(e,ref,confirmed=True)
        self.assertEqual(x['target_status'],'INCOMPLETE');self.assertEqual(x['whole_plasmid_status'],'INCOMPLETE')
    def test_full_length_with_ambiguous_backbone_cannot_pass(self):
        ref=fixture();q=list(ref.sequence);q[1200]='N';x=verify_event(event(''.join(q)),ref,confirmed=True)
        self.assertEqual(x['whole_plasmid_status'],'INCOMPLETE')

class UnassignedReportingTests(unittest.TestCase):
    def test_unassigned_events_remain_visible_in_report(self):
        from verification_align import unassigned_result
        e=event('ACGT');x=unassigned_result(e)
        self.assertEqual(x['target_status'],'REVIEW')
        outputs=project_exports([x],[],[e]);self.assertIn(b'UNASSIGNED',outputs['html'])
        with zipfile.ZipFile(io.BytesIO(outputs['zip'])) as z:
            self.assertTrue(any(n.endswith('qc.json') for n in z.namelist()))
            self.assertFalse(any(n.endswith('reference.gbk') for n in z.namelist()))

class UnfilteredEvidenceTests(unittest.TestCase):
    def test_vendor_filtered_out_high_fraction_candidate_stays_review(self):
        ref=fixture();e=event(ref.sequence[482:758])
        e.variants=[{'table':'unfiltered','Pos':'100','Ref':'A','Alt':'G','Freq':'0.2','Alt_dp':'12','Filter_Status':'FILTERED','Filter_Reason':'strand bias'}]
        x=verify_event(e,ref,confirmed=True)
        self.assertEqual(x['target_status'],'REVIEW')
        self.assertEqual(x['mixture_evidence'][0]['source'],'vendor_unfiltered_candidate_requires_review')

if __name__=='__main__':unittest.main()
