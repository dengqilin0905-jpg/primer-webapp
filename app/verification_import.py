"""Bounded ZIP import; never extract uploads or interpret generated AB1 as raw quality."""
from dataclasses import dataclass, field, asdict
from hashlib import sha256
from pathlib import PurePosixPath
import io,zipfile,re,csv,struct,gzip,collections,json
from Bio import SeqIO

MAX_ARCHIVE=40_000_000
MAX_EXPANDED=150_000_000
@dataclass
class Event:
    event_id: str
    company_id: str
    sample_name: str
    processing_batch: str
    colony_suffix: str
    archive_sha256: str
    content_sha256: str
    consensus: str = ''
    consensus_sha256: str = ''
    files: list = field(default_factory=list)
    variants: list = field(default_factory=list)
    bam_metrics: dict = field(default_factory=dict)
    ab1: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)
    batch_evidence: list = field(default_factory=list)
    clone_id: str = ''
    colony_id: str = ''
    construct_identity: str = ''
    reference_sha256: str = ''

    def summary(self):
        d=asdict(self); d.pop('consensus');return d

def bam_metrics(data):
    b=gzip.GzipFile(fileobj=io.BytesIO(data)).read(MAX_EXPANDED+1) if data[:2]==b'\x1f\x8b' else data
    if len(b)>MAX_EXPANDED or b[:4]!=b'BAM\x01': raise ValueError('Invalid or excessive BAM')
    l=struct.unpack_from('<i',b,4)[0]; header=b[8:8+l].decode(errors='replace');p=8+l
    nr=struct.unpack_from('<i',b,p)[0];p+=4;refs=[]
    for _ in range(nr):
        n=struct.unpack_from('<i',b,p)[0];p+=4;name=b[p:p+n-1].decode();p+=n
        size=struct.unpack_from('<i',b,p)[0];p+=4
        if size<0 or size>1_000_000:raise ValueError('BAM reference size limit')
        refs.append((name,size))
    if len(refs)!=1:raise ValueError('Expected one vendor consensus reference')
    cov=[0]*refs[0][1];strand_cov=[[0]*refs[0][1],[0]*refs[0][1]];records=0;names=set();forward=reverse=0;dels=collections.defaultdict(lambda:[set(),set()]);read_blocks=[]
    while p<len(b):
        n=struct.unpack_from('<i',b,p)[0];p+=4
        if n<32 or p+n>len(b):raise ValueError('Truncated BAM record')
        v=b[p:p+n];p+=n;records+=1
        rid,rp,bmq,fnc,ls,_,_,_=struct.unpack_from('<iiIIiiii',v,0)
        flag=fnc>>16;nc=fnc&65535;nl=bmq&255
        if flag&(4|256|2048):continue
        if rid!=0:continue
        read=v[32:32+nl-1].decode(errors='replace');names.add(read)
        strand=bool(flag&16);reverse+=strand;forward+=not strand
        seq_offset=32+nl+nc*4
        packed=v[seq_offset:seq_offset+(ls+1)//2]
        alphabet='=ACMGRSVTWYHKDBN'
        read_sequence=''.join(alphabet[x>>4]+alphabet[x&15] for x in packed)[:ls]
        query_pos=0
        for j in range(nc):
            cg=struct.unpack_from('<I',v,32+nl+j*4)[0];ln,op=cg>>4,cg&15
            if op in (0,7,8):
                read_blocks.append({'start0':rp,'end0':rp+ln,'sequence':read_sequence[query_pos:query_pos+ln],
                                    'read_name':read,'reverse':strand})
                for k in range(max(0,rp),min(len(cov),rp+ln)):
                    cov[k]+=1;strand_cov[strand][k]+=1
            elif op==2:dels[(rp,ln)][strand].add(read)
            if op in (0,2,3,7,8):rp+=ln
            if op in (0,1,4,7,8):query_pos+=ln
    return {'reference_name':refs[0][0],'reference_length':len(cov),'header':header,
            'self_consensus_alignment':'.consensus.insert.fasta' in header,
            'records':records,'unique_read_names':len(names),'forward':forward,'reverse':reverse,
            'read_blocks':read_blocks,'depth':cov,'depth_forward':strand_cov[0],'depth_reverse':strand_cov[1],'large_deletions':[{'start0':x,'length':n,'forward_reads':len(v[0]),
            'reverse_reads':len(v[1]),'unique_reads':len(v[0]|v[1]),
            'fraction':len(v[0]|v[1])/max(1,cov[x] if 0<=x<len(cov) else 0)} for (x,n),v in dels.items()]}

def parse_zip(data):
    if len(data)>MAX_ARCHIVE:raise ValueError('ZIP exceeds 40 MB limit')
    zh=sha256(data).hexdigest()
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        members=[x for x in z.infolist() if not x.is_dir()]
        if len(members)>500 or sum(x.file_size for x in members)>MAX_EXPANDED:raise ValueError('ZIP expansion limit')
        content=[]; payload={}
        for info in members:
            p=PurePosixPath(info.filename)
            if p.is_absolute() or '..' in p.parts or '\\' in info.filename:raise ValueError('Unsafe ZIP member')
            if info.file_size>MAX_ARCHIVE or info.file_size/max(info.compress_size,1)>1000:raise ValueError('ZIP member expansion limit')
            b=z.read(info);payload[info.filename]=b;content.append((PurePosixPath(info.filename).name,sha256(b).hexdigest()))
        if len(payload)!=len(members):raise ValueError('Duplicate ZIP member names')
        evidence=sha256(json.dumps(sorted(content)).encode()).hexdigest()
        groups=collections.defaultdict(list)
        for name,b in payload.items():
            m=re.search(r'(S\d+)-(.+?)(?:\.consensus\.insert\.fasta|_\d+\.ab1|\.filt\.var\.xls|\.var\.xls|\.bam(?:\.bai)?|\.base\.png|\.coverage\.png|\.pla_peak_check\.png|_plasmidviewer\.png)$',PurePosixPath(name).name)
            if m:groups[(m[1],m[2])].append((name,b))
        if not groups:raise ValueError('No recognized CWBIO sample records')
        events=[]
        for (company,sample),items in groups.items():
            eh=sha256(json.dumps(sorted((PurePosixPath(n).name,sha256(b).hexdigest()) for n,b in items)).encode()).hexdigest()
            e=Event(eh,company,sample,'unknown',sample.rsplit('-',1)[-1] if re.search(r'-\d+$',sample) else '',zh,evidence)
            for name,b in items:
                e.files.append({'name':name,'bytes':len(b),'sha256':sha256(b).hexdigest()})
                if name.endswith('.fasta'):
                    recs=list(SeqIO.parse(io.StringIO(b.decode('utf-8-sig')),'fasta'))
                    if len(recs)!=1:raise ValueError('Exactly one consensus per sample required')
                    q=str(recs[0].seq).upper()
                    if set(q)-set('ACGTRYSWKMBDHVN') or len(q)>50000:raise ValueError('Invalid consensus')
                    if e.consensus and e.consensus!=q:raise ValueError('Conflicting consensus records for same company sample')
                    e.consensus=q;e.consensus_sha256=sha256(q.encode()).hexdigest()
                elif name.endswith('.ab1'):
                    rec=SeqIO.read(io.BytesIO(b),'abi');raw=rec.annotations.get('abif_raw',{});comment=str(raw.get('CMNT1',''))
                    qs=rec.letter_annotations.get('phred_quality',[]);generated='Creator' in comment or 'Generated' in comment
                    e.ab1={'generated':generated,'independent_quality':not generated,'constant_q50':bool(qs) and set(qs)=={50},'comment':comment,'sequence_sha256':sha256(str(rec.seq).upper().encode()).hexdigest()}
                    e.batch_evidence.append({'source':'AB1_RUND','value':str(raw.get('RUND1','unknown'))})
                    if e.processing_batch=='unknown':e.processing_batch=str(raw.get('RUND1','unknown'))
                elif name.endswith('.bam'):
                    e.bam_metrics=bam_metrics(b)
                    m=re.search(r'/((?:20\d{6})-[^/]+)/',e.bam_metrics['header'])
                    if m:
                        e.processing_batch=m[1];e.batch_evidence.append({'source':'BAM_PG_path','value':m[1]})
                elif name.endswith('.xls'):
                    try:t=b.decode('utf-8-sig')
                    except UnicodeDecodeError:t=b.decode('gb18030')
                    for row in csv.DictReader(io.StringIO(t),delimiter='\t'):
                        e.variants.append({**row,'table':'filtered' if '.filt.' in name else 'unfiltered'})
            if not e.consensus:e.warnings.append('NO_SEQUENCE_IMAGES_ONLY')
            if e.ab1.get('generated'):e.warnings.append('GENERATED_AB1_NOT_INDEPENDENT_QUALITY')
            if e.ab1 and e.consensus_sha256!=e.ab1['sequence_sha256']:e.warnings.append('AB1_CONSENSUS_CONFLICT')
            if e.processing_batch=='unknown':
                dates=[x.date_time for x in members if x.filename in {n for n,_ in items}]
                if dates:e.processing_batch='zip-metadata:'+ '-'.join(map(str,min(dates)[:3]))
            events.append(e)
        return events

def merge_events(existing, incoming):
    """Event identity includes evidence hash; sample names never act as storage keys."""
    out={e.event_id:e for e in existing}
    for e in incoming:out.setdefault(e.event_id,e)
    grouped=collections.defaultdict(list)
    for e in out.values():grouped[e.sample_name].append(e)
    for values in grouped.values():
        if len(values)>1:
            for e in values:
                if 'SAMPLE_NAME_COLLISION' not in e.warnings:e.warnings.append('SAMPLE_NAME_COLLISION')
    return list(out.values())
