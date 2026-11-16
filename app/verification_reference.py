"""Immutable, assembly-backed reference snapshots. Coordinates are 0-based half-open."""
from dataclasses import dataclass, asdict
from hashlib import sha256
import io, json
from Bio import SeqIO
from Bio.Seq import Seq

@dataclass(frozen=True)
class Region:
    name: str
    start: int
    end: int
    kind: str
    required: bool = True

@dataclass(frozen=True)
class Reference:
    construct: str
    operational_label: str
    sequence: str
    sequence_sha256: str
    genbank: str
    regions: tuple[Region, ...]
    provenance_json: str
    aliases: tuple[str, ...] = ()
    gene_id: str = ''
    canonical_identity: str = ''

    def snapshot(self):
        return asdict(self)

def freeze_reference(record, *, backbone_sequence, replace_start, replace_end,
                     insert_sequence, source, operational_label, aliases=(),
                     gene_id='', canonical_identity='', required_features=()):
    """Never infer insert extent from a CDS feature. Verify the complete assembly."""
    v=str(backbone_sequence).upper(); ins=str(insert_sequence).upper(); s=str(record.seq).upper()
    a,b=int(replace_start),int(replace_end)
    if not source or not 0 <= a <= b <= len(v) or not ins:
        raise ValueError('Explicit assembly provenance and nonempty insert are required')
    if set(v+ins)-set('ACGT') or v[:a]+ins+v[b:] != s:
        raise ValueError('Complete reference does not equal provenance-backed assembly')
    end=a+len(ins)
    if a<18 or end+18>len(s):
        raise ValueError('Junction flanks need an explicit circular boundary model')
    regions=[Region('design_insert',a,end,'insert'),Region('junction_5',a-18,min(a+18,end),'junction'),
             Region('junction_3',max(a,end-18),end+18,'junction')]
    tags=[]
    for f in record.features:
        label=';'.join(f.qualifiers.get('label',f.qualifiers.get('gene',[])))
        for part in f.location.parts:
            x,y=int(part.start),int(part.end)
            is_tag=label.lower().replace(' ','') in ('3xha','ha','flag','gfp','nyfp','cyfp')
            if is_tag: tags.append((x,y,label))
            required=(f.type=='CDS' and x>=a and y<=end) or is_tag or label in required_features
            if required: regions.append(Region(label or f.type,x,y,'tag' if is_tag else f.type))
            if f.type=='CDS' and x>=a and y<=end and y-x>=3 and (y-x)%3==0 and s[y-3:y] in ('TAA','TAG','TGA'):
                regions.append(Region((label or 'CDS')+'_terminal_stop',y-3,y,'stop'))
    for x,y,label in tags:
        if x>=end and (x-end)%3==0:
            if x>end: regions.append(Region('fusion_linker',end,x,'linker'))
            stop=next((i for i in range(y,len(s)-2,3) if s[i:i+3] in ('TAA','TAG','TGA')),None)
            if stop is None: raise ValueError('Required fusion stop cannot be recovered')
            regions.append(Region('fusion_stop',y,stop+3,'stop'))
    text=io.StringIO(); SeqIO.write(record,text,'genbank')
    provenance={'source':source,'backbone_sha256':sha256(v.encode()).hexdigest(),
                'original_replacement':[a,b],'insert_interval':[a,end],
                'insert_sha256':sha256(ins.encode()).hexdigest(),'assembly_verified':True,
                'junction_5':s[a-18:a+18],'junction_3':s[end-18:end+18]}
    return Reference(record.id,operational_label,s,sha256(s.encode()).hexdigest(),text.getvalue(),
                     tuple(regions),json.dumps(provenance,sort_keys=True),tuple(aliases),gene_id,canonical_identity)

def recover_reference(record, original, *, source, operational_label, aliases=(), **kwargs):
    """Restore from mapped arms plus original vector, then validate every retained base."""
    fs=[f for f in record.features if f.qualifiers.get('label')==['F_homology_arm']]
    rs=[f for f in record.features if f.qualifiers.get('label')==['R_homology_arm']]
    if len(fs)!=1 or len(rs)!=1: raise ValueError('Unique mapped F/R arms required')
    f,rr=fs[0],rs[0]; s=str(record.seq).upper(); v=str(original.seq).upper()
    a=int(f.location.end); e=int(rr.location.start)
    if a>e: raise ValueError('Origin-crossing insert requires explicit assembly provenance')
    left=s[int(f.location.start):a]; right=s[e:int(rr.location.end)]
    def unique(needle):
        hits=[i for i in range(len(v)) if v.startswith(needle,i)]
        if len(hits)!=1: raise ValueError('Arms are not unique in original backbone')
        return hits[0]
    old_a=unique(left)+len(left); old_b=unique(right)
    if old_a!=a: raise ValueError('Reference origin differs from assembly provenance')
    return freeze_reference(record,backbone_sequence=v,replace_start=old_a,replace_end=old_b,
                            insert_sequence=s[a:e],source=source,operational_label=operational_label,
                            aliases=aliases,**kwargs)


def require_features(reference, labels):
    """Freeze additional mandatory map elements without changing reference sequence."""
    from dataclasses import replace
    record=SeqIO.read(io.StringIO(reference.genbank),'genbank')
    regions=list(reference.regions);found=set()
    for f in record.features:
        label=';'.join(f.qualifiers.get('label',f.qualifiers.get('gene',[])))
        if label in labels:
            found.add(label)
            for part in f.location.parts:
                reg=Region(label,int(part.start),int(part.end),f.type)
                if reg not in regions:regions.append(reg)
    if set(labels)-found:raise ValueError('Required feature not present in reference map')
    provenance=json.loads(reference.provenance_json);provenance['additional_required_features']=sorted(labels)
    return replace(reference,regions=tuple(regions),provenance_json=json.dumps(provenance,sort_keys=True))


def reference_from_canonical(item):
    """Consume canonical service inputs; never redesign or reconstruct business naming."""
    from genbank_export import build_genbank_record, _resolve_cut_window
    row=dict(item.map_input);v=str(row['vector_sequence']).upper()
    a,b,_=_resolve_cut_window(row,v)
    return freeze_reference(build_genbank_record(row),backbone_sequence=v,replace_start=a,replace_end=b,
        insert_sequence=row['insert_seq'],source=json.dumps(asdict(item.provenance),sort_keys=True),
        operational_label=item.operational_label,gene_id=item.gene_id,
        canonical_identity=item.canonical_identity)
