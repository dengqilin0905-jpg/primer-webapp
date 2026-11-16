"""Read-only project-map catalog; assembly equality is mandatory for every reference."""
from pathlib import Path
from resource_config import resource_path
from identity_audit import RICE_IDENTITY_REGISTRY
import csv,re
from Bio import SeqIO
from naming import build_operational_label
from verification_reference import recover_reference

def load_catalog(project_root,repo_root):
    root=Path(project_root);repo=Path(repo_root);vectors=[];references=[];errors=[]
    if not root.is_dir():return [],['No project directory configured; upload a reference and explicit assembly provenance.']
    vector_root=resource_path('backbone_root',repo)
    for path in list((vector_root/'resources/backbones').glob('*.dna'))+[p for d in root.iterdir() if d.is_dir() for p in d.glob('*.dna')]:
        try:vectors.append((path,SeqIO.read(path,'snapgene')))
        except Exception as exc:errors.append(f'{path.name}: {type(exc).__name__}')
    registry=RICE_IDENTITY_REGISTRY;identities={}
    if registry.is_file():
        for row in csv.DictReader(registry.open(),delimiter='\t'):
            identities.setdefault(row.get('gene_ID',''),[]).append(row.get('construct_ID',''))
    paths=[]
    for d in root.iterdir():
        if not d.is_dir():continue
        for rel in ('results/maps','primer_output/construct_maps','rice_chip_2/plasmids'):
            folder=d/rel
            if folder.is_dir():paths.extend(p for p in folder.iterdir() if p.suffix.lower() in ('.gb','.gbk'))
    seen=set()
    for path in sorted(paths):
        try:
            record=SeqIO.read(path,'genbank');name=record.id
            fs=[f for f in record.features if f.qualifiers.get('label')==['F_homology_arm']]
            rs=[f for f in record.features if f.qualifiers.get('label')==['R_homology_arm']]
            if len(fs)!=1 or len(rs)!=1:continue
            a,e=int(fs[0].location.end),int(rs[0].location.start)
            targets=[f for f in record.features if f.type=='CDS' and int(f.location.start)>=a and int(f.location.end)<=e and 'locus_tag' in f.qualifiers]
            gene=targets[0].qualifiers['locus_tag'][0] if len(targets)==1 else ''
            backbone=name.split('_',1)[0]
            mid=re.search(r'(?:^|_)(?:Os)?(RTF\d+|TF\d+)(?:_|$)',name);identity=mid[1] if mid else ''
            if identity and gene and identity not in identities.get(gene,[]):
                errors.append(f'{name}: canonical identity registry conflict');continue
            symbol=(targets[0].qualifiers.get('gene') or targets[0].qualifiers.get('label'))[0] if targets else name.split('_',1)[-1]
            if identity:symbol=re.sub(r'^(?:RTF|TF)\d+_','',symbol)
            species='Os' if gene.startswith('LOC_Os') else None
            try:label=build_operational_label(symbol,species=species,identifier=identity or None,backbone=backbone)
            except ValueError:
                errors.append(f'{name}: explicit species/identity metadata required');continue
            for vp,vec in vectors:
                try:
                    ref=recover_reference(record,vec,source=f'{path}; original={vp}',operational_label=label,
                                          aliases=(name,),gene_id=gene,canonical_identity=identity)
                    key=(ref.construct,ref.sequence_sha256)
                    if key not in seen:references.append(ref);seen.add(key)
                    break
                except ValueError:continue
            else:errors.append(f'{name}: no exact original-backbone assembly provenance')
        except Exception as exc:errors.append(f'{path.name}: {exc}')
    return references,errors
