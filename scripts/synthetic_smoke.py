"""Self-contained synthetic Contract/Web/Batch/export acceptance; no lab resources."""
import io,json,os,random,sys,tempfile,hashlib,zipfile
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
ROOT=Path(__file__).resolve().parents[1]

def prepare(folder):
    folder=Path(folder);(folder/'config').mkdir(parents=True);(folder/'resources').mkdir()
    rng=random.Random(904);v=''.join(rng.choice('ACGT') for _ in range(1800));a,b=400,460
    v=v[:a-6]+'ACTAGT'+v[a:b]+'CCCGGG'+v[b+6:]
    v=v[:b+18]+'TAC'*30+'TAA'+v[b+111:]
    rec=SeqRecord(Seq(v),id='SYNTHETIC_VECTOR',description='Artificial software fixture; no experimental provenance')
    rec.annotations={'molecule_type':'DNA','topology':'circular'}
    rec.features=[SeqFeature(FeatureLocation(a,b),type='misc_feature',qualifiers={'label':['synthetic_placeholder']}),SeqFeature(FeatureLocation(b+18,b+108),type='misc_feature',qualifiers={'label':['synthetic_tag']})]
    SeqIO.write(rec,folder/'resources/vector.gbk','genbank')
    hf=v[a-18:a];hr=str(Seq(v[b:b+18]).reverse_complement())
    preset={'preset_id':'SYN','preset_version':1,'display_name':'Synthetic demonstration only','aliases':[],
      'status':'current','validation_level':'computational','approved_by_user':True,
      'vector':{'path':'resources/vector.gbk','sha256':hashlib.sha256(v.encode()).hexdigest(),'sequence_length':len(v),'topology':'circular','source_format':'genbank'},
      'cloning':{'insertion_mode':'restriction_site_join','restriction_site_f':{'enzyme':'SpeI','sequence':'ACTAGT'},'restriction_site_r':{'enzyme':'SmaI','sequence':'CCCGGG'},
        'homology':{'arm_length':18,'homo_f':hf,'homo_r':hr},'replacement':{'source_feature':'synthetic_placeholder','boundary_start':a,'boundary_end':b},
        'fusion':{'order':'target -> synthetic_tag','architecture':'C-terminal synthetic tag','construct_type':'C_terminal_tag'}},
      'design_defaults':{'remove_stop':True},'qc':{'primer_engine_frame_check':False,'fusion_frame_validation':True},'provenance':{'source_project':'synthetic demo','validation_note':'Artificial software fixture, never a laboratory preset'}}
    (folder/'config/presets.json').write_text(json.dumps({'registry_version':1,'current_versions':{'SYN':1},'presets':[preset]}))
    (folder/'identity.tsv').write_text('gene_id\tcanonical_identity\tphysical_plasmid_status\nLOC_SYNTHETIC\tTF999\tmissing\n')
    (folder/'annotation.tsv').write_text('geneID\tgene_symbol\tAtBestHit\tArabidopsis_symbol\nLOC_SYNTHETIC\tOsSynthetic\t\t\n')
    (folder/'curated.tsv').write_text('Gene_ID\tGene_symbol\n')
    config={'backbone_registry':str(folder/'config/presets.json'),'backbone_root':str(folder),'identity_registry':str(folder/'identity.tsv'),
            'rice_annotation':str(folder/'annotation.tsv'),'curated_annotation':str(folder/'curated.tsv'),'projects_root':str(folder/'projects')}
    (folder/'private-config.json').write_text(json.dumps(config))
    return config,hf

def run(folder):
    folder=Path(folder);config,hf=prepare(folder)
    os.environ['PRIMER_PRIVATE_CONFIG']=str(folder/'private-config.json')
    sys.path.insert(0,str(ROOT/'app'))
    from project_contract import CANONICAL_OUTPUTS
    from project_orchestration import run_project_directory,build_project_package,build_contract_artifacts
    from project_validator import audit_project
    from verification_reference import reference_from_canonical
    from verification_import import Event
    from verification_align import verify_event
    from verification_export import project_exports
    project=folder/'projects/synthetic_contract';(project/'inputs').mkdir(parents=True)
    manifest={'contract_version':1,'project_id':project.name,'project_type':'primer_design','created_at':'2026-01-01T00:00:00Z','species':['Rice'],
      'constructs':[{'gene_id':'LOC_SYNTHETIC','symbol':'OsSynthetic','canonical_identity':'TF999','operational_label':'SYN_OsTF999','full_construct_name':'SYN_TF999_OsSynthetic_synthetic_tag',
       'backbone':'SYN','preset_id':'SYN','preset_version':1,'remove_stop':True,'physical_plasmid_status':'missing','overrides':[],'status':'planned'}],
      'outputs':dict(CANONICAL_OUTPUTS),'provenance':{'identity_sources':[str(folder/'identity.tsv')],'vector_hash':json.loads((folder/'config/presets.json').read_text())['presets'][0]['vector']['sha256'],'template_source':'generic software workbook'}}
    runtime={'sequences':{'LOC_SYNTHETIC':'ATG'+'GCC'*100+'TAA'},'identity_sources':{'identity_registry':[str(folder/'identity.tsv')]},
      'structural_qc':{'LOC_SYNTHETIC':{kind:{'status':'PASS','detail':'synthetic test evidence'} for kind in ['hairpin','homodimer','heterodimer']}},
      'sequencing_inputs':{'LOC_SYNTHETIC':{'reference_primer_direction':'F','reference_sequence':hf,'resistance':'AmpR'}}}
    # Use an exact template-legal value in the canonical control input; alias handling has independent tests.
    runtime['sequencing_inputs']['LOC_SYNTHETIC']['resistance']='氨苄'
    (project/'project.yaml').write_text(json.dumps(manifest));(project/'inputs/runtime_inputs.json').write_text(json.dumps(runtime));(project/'WORK_LOG.md').write_text('# Synthetic fixture\n')
    batch=run_project_directory(project)
    package,web=build_project_package((project/'project.yaml').read_bytes(),(project/'inputs/runtime_inputs.json').read_bytes(),(project/'WORK_LOG.md').read_bytes(),active_projects_root=folder/'projects',exclude_project_root=project)
    expected=build_contract_artifacts(manifest,runtime,active_projects_root=folder/'projects',exclude_project_root=project)
    assert batch['conformance']['status']==web['conformance']['status']=='PASS'
    assert batch['normalized_canonical_result']==web['normalized_canonical_result']==expected['normalized_canonical_result']
    assert audit_project(project)['status']=='PASS'
    with zipfile.ZipFile(io.BytesIO(package)) as z:assert 'results/sequencing/sequencing_order.xlsx' in z.namelist()
    # Independent reference/consensus verification uses the released canonical assembly.
    item=expected['canonical_result'].constructs[0]
    ref=reference_from_canonical(item)
    q=ref.sequence;h=hashlib.sha256(q.encode()).hexdigest()
    event=Event(h,'S_SYNTHETIC','SYN_OsTF999-1','synthetic','1',h,h,consensus=q,consensus_sha256=h,
       bam_metrics={'reference_length':len(q),'depth':[20]*len(q),'depth_forward':[10]*len(q),'depth_reverse':[10]*len(q),'forward':10,'reverse':10})
    result=verify_event(event,ref,confirmed=True);assert result['target_status']=='PASS'
    outputs=project_exports([result],[ref],[event]);assert b'Actual aligned bases' in outputs['html']
    assert result['differences']==[] and all(x['coverage_fraction']==1 for x in result['regions'])
    print(json.dumps({'synthetic_contract':'PASS','web_batch_parity':'PASS','target':'PASS','required_coverage':1.0,'native_snapgene_alignment':False}))
if __name__=='__main__':
    with tempfile.TemporaryDirectory(prefix='primer-synthetic-') as folder:run(folder)
