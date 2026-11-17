"""Portable distribution, public configuration and synthetic workflow acceptance."""
import io,json,os,subprocess,sys,tempfile,unittest,hashlib,struct
from pathlib import Path
from unittest.mock import patch
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
from resource_config import resource_path,private_settings,REPO_ROOT
from backbone_presets import list_backbone_presets,get_backbone_preset,BackbonePresetError
from project_contract import minimal_backbone_name
from verification_catalog import load_catalog
from seq_io import parse_uploaded_records

class PrivateConfigurationTests(unittest.TestCase):
    def test_default_has_no_lab_presets(self):self.assertEqual(list_backbone_presets(),[])
    def test_missing_project_directory_does_not_crash_catalog(self):
        with tempfile.TemporaryDirectory() as d:
            refs,issues=load_catalog(Path(d)/'missing',REPO_ROOT)
            self.assertEqual(refs,[]);self.assertTrue(issues)
    def test_relative_paths_are_relative_to_external_config(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'resources.json';p.write_text(json.dumps({'projects_root':'projects'}))
            with patch.dict(os.environ,{'PRIMER_PRIVATE_CONFIG':str(p)}):self.assertEqual(resource_path('projects_root','fallback'),Path(d)/'projects')
    def test_environment_overrides_external_file(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'resources.json';p.write_text(json.dumps({'projects_root':'first'}))
            with patch.dict(os.environ,{'PRIMER_PRIVATE_CONFIG':str(p),'PRIMER_PROJECTS_ROOT':str(Path(d)/'second')}):self.assertEqual(resource_path('projects_root','fallback'),Path(d)/'second')
    def test_invalid_config_does_not_silently_use_another_project(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'resources.json';p.write_text('[]')
            with patch.dict(os.environ,{'PRIMER_PRIVATE_CONFIG':str(p)}):
                with self.assertRaises(ValueError):private_settings()
    def test_simple_custom_backbone_is_an_explicit_name(self):
        self.assertEqual(minimal_backbone_name('SYN'),'SYN')
        with self.assertRaises(ValueError):minimal_backbone_name('../../vector')
    def test_external_genbank_preset_and_hash_validation(self):
        sys.path.insert(0,str(REPO_ROOT/'scripts'));from synthetic_smoke import prepare
        with tempfile.TemporaryDirectory() as d:
            prepare(d);p=Path(d)/'config/presets.json';before=p.read_bytes()
            preset=get_backbone_preset('SYN',path=p);self.assertEqual(preset['cloning']['homology']['arm_length'],18)
            from backbone_presets import validate_backbone_preset
            preset['vector']['sha256']='0'*64
            with self.assertRaises(BackbonePresetError):validate_backbone_preset(preset,registry_path=p)
            self.assertEqual(before,p.read_bytes())
    def test_genbank_upload_preserves_features_and_topology(self):
        rec=SeqRecord(Seq('ATG'+'GCT'*20+'TAA'),id='SYNTHETIC',description='artificial fixture')
        rec.annotations={'molecule_type':'DNA','topology':'circular'};rec.features=[SeqFeature(FeatureLocation(0,len(rec)),type='CDS',qualifiers={'label':['synthetic_CDS']})]
        b=io.StringIO();SeqIO.write(rec,b,'genbank');upload=io.BytesIO(b.getvalue().encode());upload.name='synthetic.gbk'
        parsed=parse_uploaded_records(upload)[0];self.assertEqual(parsed['sequence'],str(rec.seq));self.assertEqual(parsed['topology'],'circular');self.assertTrue(parsed['features'])
    def test_synthetic_snapgene_stays_binary_and_uses_primary_parser(self):
        q='ACGT'*30
        header=struct.pack('>BI8sHHH',9,14,b'SnapGene',1,1,1)
        dna=struct.pack('>BI',0,len(q)+1)+b'\x01'+q.encode()
        upload=io.BytesIO(header+dna);upload.name='synthetic.dna'
        with patch('snapgene_reader.snapgene_file_to_dict',side_effect=AssertionError('Primary binary parser must succeed')):
            parsed=parse_uploaded_records(upload)[0]
        self.assertEqual(parsed['sequence'],q);self.assertEqual(parsed['topology'],'circular')
    def test_external_asset_root_is_the_same_for_validation_and_loading(self):
        sys.path.insert(0,str(REPO_ROOT/'scripts'));from synthetic_smoke import prepare
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);prepare(d);extra=d/'separate-config';extra.mkdir()
            registry=extra/'registry.json';registry.write_bytes((d/'config/presets.json').read_bytes())
            cfg=d/'override.json';cfg.write_text(json.dumps({'backbone_registry':str(registry),'backbone_root':str(d)}))
            script="from backbone_presets import get_backbone_preset;from design_config import shared_preset_design_config;import hashlib;p=get_backbone_preset('SYN');c=shared_preset_design_config('SYN');assert hashlib.sha256(c['vector_sequence'].encode()).hexdigest()==p['vector']['sha256']"
            env=dict(os.environ,PRIMER_PRIVATE_CONFIG=str(cfg),PYTHONPATH=str(REPO_ROOT/'app'))
            result=subprocess.run([sys.executable,'-B','-c',script],env=env,capture_output=True,text=True,timeout=20)
            self.assertEqual(result.returncode,0,result.stderr)
    def test_synthetic_contract_web_batch_and_verification(self):
        env=dict(os.environ);env.pop('PRIMER_PRIVATE_CONFIG',None)
        p=subprocess.run([sys.executable,'-B',str(REPO_ROOT/'scripts/synthetic_smoke.py')],text=True,capture_output=True,env=env,timeout=40)
        self.assertEqual(p.returncode,0,p.stderr);self.assertIn('"target": "PASS"',p.stdout);self.assertIn('"web_batch_parity": "PASS"',p.stdout)

class StandaloneAppTests(unittest.TestCase):
    def test_registration_login_all_tabs_and_custom_backbone_without_presets(self):
        from streamlit.testing.v1 import AppTest
        from db import get_conn
        at=AppTest.from_file(str(REPO_ROOT/'app/app.py'),default_timeout=30).run();self.assertFalse(at.exception)
        at.text_input(key='reg_u').input('synthetic_user');at.text_input(key='reg_p1').input('synthetic-password');at.text_input(key='reg_p2').input('synthetic-password')
        next(b for b in at.button if b.label=='Create account').click().run();self.assertFalse(at.exception)
        at.text_input(key='login_u').input('synthetic_user');at.text_input(key='login_p').input('synthetic-password');next(b for b in at.button if b.label=='Login').click().run()
        self.assertFalse(at.exception);self.assertTrue(any(t.label=='Sequencing Verification' for t in at.tabs))
        uid=at.session_state['user']['id'];conn=get_conn()
        record=SeqIO.read(REPO_ROOT/'examples/synthetic_vector.gbk','genbank');vector=str(record.seq)
        from seq_io import _feature_to_json
        conn.execute("INSERT INTO backbones(user_id,name,enzyme_sites,homo_f,homo_r,vector_sequence,insertion_mode,restriction_site_f,restriction_site_r,vector_features_json,vector_topology,created_at) VALUES(?,?,?,?,?,?,?,?,?,?,?,?)",(uid,'SYN','SpeI/SmaI',vector[382:400],str(Seq(vector[460:478]).reverse_complement()),vector,'restriction_site_join','ACTAGT','CCCGGG',json.dumps([_feature_to_json(f) for f in record.features]),'circular','synthetic'))
        conn.execute('INSERT INTO genes(user_id,gene_name,sequence,created_at) VALUES(?,?,?,?)',(uid,'OsSynthetic','ATG'+'GCT'*99+'TAA','synthetic'));conn.commit();conn.close()
        at.run();self.assertFalse(at.exception);self.assertEqual(at.radio(key='design_backbone_source').value,'My/custom backbone')
        at.multiselect(key='design_selected_genes').set_value(['OsSynthetic']).run()
        at.button(key='design_add_batch').click().run();self.assertFalse(at.exception)
        at.button(key='design_run').click().run();self.assertFalse(at.exception)
        self.assertTrue(at.session_state['last_run_id']);self.assertFalse(at.session_state['last_failures'])
        self.assertGreater(len(at.get('download_button')),0)
    def test_ci_cannot_deploy_on_main_push(self):
        self.assertFalse((REPO_ROOT/'.github/workflows/deploy.yml').exists())
        s=(REPO_ROOT/'.github/workflows/software-ci.yml').read_text();self.assertNotIn('push:',s);self.assertNotIn('secrets.',s);self.assertNotIn('remote_deploy',s)
        s=(REPO_ROOT/'scripts/remote_deploy.sh').read_text();self.assertNotIn('rm -rf',s);self.assertNotIn('systemctl',s)

if __name__=='__main__':unittest.main()
