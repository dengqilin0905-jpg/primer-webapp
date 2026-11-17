"""Pre-release acceptance tests; all sequence content here is synthetic and public."""
import unittest,random,io,json,zipfile,hashlib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
from openpyxl import load_workbook
from verification_reference import freeze_reference
from verification_align import verify_event
from verification_import import Event
from verification_export import project_exports

def fusion_fixture():
    rng=random.Random(20261116)
    v=''.join(rng.choice('ACGT') for _ in range(1800));a,b=400,460
    # Insert comprises 5-prime non-CDS, CDS and extra insert sequence.
    prefix='CGTACGTACGTA';cds='ATG'+'GCT'*39;tail='GGTGGT';insert=prefix+cds+tail
    linker='GGT'*6;tag='TAC'*30
    v=v[:b]+linker+tag+'TAA'+v[b+len(linker)+len(tag)+3:]
    sequence=v[:a]+insert+v[b:];end=a+len(insert);tag_start=end+len(linker);tag_end=tag_start+len(tag)
    rec=SeqRecord(Seq(sequence),id='SYNTHETIC_FUSION',description='public synthetic acceptance fixture')
    rec.annotations={'molecule_type':'DNA','topology':'circular'}
    rec.features=[SeqFeature(FeatureLocation(a+len(prefix),a+len(prefix)+len(cds)),type='CDS',qualifiers={'label':['synthetic_CDS']}),
        SeqFeature(FeatureLocation(a,a+len(prefix)),type='misc_feature',qualifiers={'label':['synthetic_non_CDS']}),
        SeqFeature(FeatureLocation(tag_start,tag_end),type='CDS',qualifiers={'label':['3xHA']})]
    ref=freeze_reference(rec,backbone_sequence=v,replace_start=a,replace_end=b,insert_sequence=insert,
        source='synthetic assembly fixture',operational_label='SYNTHETIC_OsTF999')
    return ref,a-18,tag_end+3

def supported_event(q):
    h=hashlib.sha256(q.encode()).hexdigest()
    return Event(h,'S_SYNTHETIC','SYNTHETIC_OsTF999-1','synthetic','1',h,h,
        consensus=q,consensus_sha256=h,bam_metrics={'reference_length':len(q),'depth':[20]*len(q),
        'depth_forward':[10]*len(q),'depth_reverse':[10]*len(q),'forward':10,'reverse':10})

class PreReleaseAcceptanceTests(unittest.TestCase):
    def test_complete_insert_non_cds_junction_tag_stop_pass_and_exports(self):
        ref,start,end=fusion_fixture();e=supported_event(ref.sequence[start:end]);x=verify_event(e,ref,confirmed=True)
        self.assertEqual(x['target_status'],'PASS');self.assertEqual(x['whole_plasmid_status'],'INCOMPLETE')
        self.assertEqual(x['differences'],[])
        self.assertTrue({'insert','junction','CDS','tag','linker','stop'} <= {row['kind'] for row in x['regions']})
        self.assertTrue(all(row['coverage_fraction']==1.0 and row['status']=='PASS' for row in x['regions']))
        out=project_exports([x],[ref],[e]);wb=load_workbook(io.BytesIO(out['xlsx']))
        headers=[c.value for c in wb['Samples'][1]]
        self.assertEqual(wb['Samples'].cell(2,headers.index('target_status')+1).value,'PASS')
        self.assertEqual(wb['Differences'].max_row,1)
        with zipfile.ZipFile(io.BytesIO(out['zip'])) as z:
            gb=next(n for n in z.namelist() if n.endswith('reference.gbk'))
            exported=SeqIO.read(io.StringIO(z.read(gb).decode()),'genbank')
            self.assertEqual(hashlib.sha256(str(exported.seq).upper().encode()).hexdigest(),ref.sequence_sha256)
            self.assertTrue(any(f.qualifiers.get('label')==['3xHA'] for f in exported.features))
            snap=json.loads(z.read(next(n for n in z.namelist() if n.endswith('reference_snapshot.json'))))
            self.assertEqual(snap['sequence_sha256'],x['reference_sha256'])
            self.assertEqual(hashlib.sha256(json.dumps(snap,sort_keys=True).encode()).hexdigest(),x['reference_snapshot_sha256'])
            self.assertEqual(z.read(next(n for n in z.namelist() if n.endswith('differences.tsv'))),b'')
    def test_non_cds_insert_change_prevents_pass(self):
        ref,start,end=fusion_fixture();q=list(ref.sequence[start:end]);q[20]='A' if q[20]!='A' else 'C'
        self.assertEqual(verify_event(supported_event(''.join(q)),ref,confirmed=True)['target_status'],'FAIL')
    def test_tag_mutation_prevents_pass(self):
        ref,start,end=fusion_fixture();tag=next(row for row in ref.regions if row.kind=='tag')
        q=list(ref.sequence[start:end]);j=tag.start-start+10;q[j]='G' if q[j]!='G' else 'A'
        self.assertEqual(verify_event(supported_event(''.join(q)),ref,confirmed=True)['target_status'],'FAIL')
    def test_cds_without_junction_or_tag_is_incomplete(self):
        ref,start,end=fusion_fixture();cds=next(row for row in ref.regions if row.kind=='CDS')
        x=verify_event(supported_event(ref.sequence[cds.start:cds.end]),ref,confirmed=True)
        self.assertEqual(x['target_status'],'INCOMPLETE')
    def test_full_fusion_reverse_orientation_passes(self):
        ref,start,end=fusion_fixture();q=str(Seq(ref.sequence[start:end]).reverse_complement())
        x=verify_event(supported_event(q),ref,confirmed=True)
        self.assertEqual(x['target_status'],'PASS');self.assertEqual(x['orientation'],'-')

if __name__=='__main__':unittest.main()
