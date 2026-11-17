"""Independent design-reference alignment and conservative regional QC."""
from collections import Counter
from hashlib import sha256
import json
from dataclasses import asdict
from Bio import Align, SeqIO
from Bio.Seq import Seq
import io

def align_consensus(reference, query):
    s=reference.sequence;n=len(s)
    if not query or n*len(query)>35_000_000:raise ValueError('Alignment size limit or empty query')
    circular=SeqIO.read(io.StringIO(reference.genbank),'genbank').annotations.get('topology')=='circular'
    best=None
    for direction,q in [('+',query),('-',str(Seq(query).reverse_complement()))]:
        offsets=Counter()
        for j in range(0,max(1,len(q)-20),25):
            k=q[j:j+21]
            if len(k)<21 or set(k)-set('ACGT'):continue
            start=0
            while True:
                i=(s+s[:20] if circular else s).find(k,start)
                if i<0:break
                offsets[(i-j)%n if circular else i-j]+=1;start=i+1
        off=offsets.most_common(1)[0][0] if offsets else 0
        margin=min(len(q),max(24,n-len(q)))
        rotation=(off-margin)%n if circular else 0
        target=(s[rotation:]+s[:rotation]+(s[rotation:rotation+48] if rotation+48<=n else (s+s)[rotation:rotation+48])) if circular else s
        a=Align.PairwiseAligner();a.mode='global';a.match_score=2;a.mismatch_score=-4;a.open_gap_score=-8;a.extend_gap_score=-.1
        # Free reference flanks; query must be independently accounted for in full.
        a.query_left_gap_score=0;a.query_right_gap_score=0
        aln=a.align(target,q)[0]
        if best is None or aln.score>best[0]:best=(aln.score,aln,q,direction,target,rotation)
    score,aln,q,direction,target,rotation=best;c=aln.coordinates
    covered={};mapped={};diffs=[];segments=[];matched=0
    for j in range(c.shape[1]-1):
        x1,x2=int(c[0,j]),int(c[0,j+1]);y1,y2=int(c[1,j]),int(c[1,j+1])
        if y1==y2 and y1 in (0,len(q)):continue
        segments.append({'reference_start0':(rotation+x1)%n,'reference_end0':(rotation+x2)%n,
                         'reference_unwrapped_start0':x1,'reference_unwrapped_end0':x2,
                         'query_start0':y1 if direction=='+' else len(q)-y2,
                         'query_end0':y2 if direction=='+' else len(q)-y1,
                         'query_oriented_start0':y1,'query_oriented_end0':y2,'orientation':direction})
        if x2>x1 and y2>y1:
            for k in range(x2-x1):
                pos=(rotation+x1+k)%n;qp=y1+k;original_qp=qp if direction=='+' else len(q)-1-qp
                covered[pos]=q[qp];mapped[pos]=original_qp
                if target[x1+k]==q[qp] and q[qp] in 'ACGT':matched+=1
                else:diffs.append({'type':'SUB','start0':pos,'length':1,'ref':target[x1+k],'alt':q[qp],'query_start0':original_qp})
        elif x2>x1:
            diffs.append({'type':'DEL','start0':(rotation+x1)%n,'length':x2-x1,'ref':target[x1:x2], 'alt':'','query_start0':y1 if direction=='+' else len(q)-y1})
        else:
            diffs.append({'type':'INS','start0':(rotation+x1)%n,'length':y2-y1,'ref':'','alt':q[y1:y2],'query_start0':y1 if direction=='+' else len(q)-y2})
    for d in diffs:
        a=d['start0'];size=d['length'] if d['type']!='INS' else 1
        names=[reg.name for reg in reference.regions if reg.required and a<reg.end and a+size>reg.start]
        d['required_regions']=names
        # Insertion at a right-hand boundary also belongs to its preceding feature.
        if d['type']=='INS':d['required_regions']+= [reg.name for reg in reference.regions if reg.required and reg.end==a]
        labels=[]
        rec=SeqIO.read(io.StringIO(reference.genbank),'genbank')
        for f in rec.features:
            if f.type=='source' or len(f)>=n:continue
            if a in f.location:labels.extend(f.qualifiers.get('label',[f.type]))
        d['features']=labels;d['classification']='required_target' if d['required_regions'] else 'retained_backbone'
    return {'orientation':direction,'score':score,'matched_bases':matched,'query_length':len(q),
            'reference_rotation0':rotation,'coverage':covered,'query_mapping':mapped,'segments':segments,'differences':diffs}

def verify_event(event,reference,*,confirmed=False,minimum_depth=5,minimum_variant_reads=3):
    result={'event_id':event.event_id,'company_id':event.company_id,'sample_name':event.sample_name,
            'batch':event.processing_batch,'clone_id':event.clone_id,'colony_id':event.colony_id,
            'colony_suffix':event.colony_suffix,'construct':reference.construct,
            'reference_sha256':reference.sequence_sha256,'consensus_sha256':event.consensus_sha256,
            'reference_snapshot_sha256':sha256(json.dumps(reference.snapshot(),sort_keys=True).encode()).hexdigest(),
            'policy_sha256':sha256(json.dumps({'regions':[asdict(x) for x in reference.regions],'min_depth':minimum_depth,'min_variant_reads':minimum_variant_reads,'mixture_fraction':.05,'algorithm_version':2,'deletion_anchor_bp':8,'actual_spanning_reads_required':True},sort_keys=True).encode()).hexdigest(),
            'policy':{'algorithm_version':2,'deletion_anchor_bp':8,'actual_spanning_reads_required':True,'minimum_depth':minimum_depth,'minimum_variant_reads':minimum_variant_reads,'mixture_fraction':.05,'junction_flank_bp':18,'per_base_bidirectional_support':True},
            'target_status':'INCOMPLETE','whole_plasmid_status':'INCOMPLETE','regions':[],
            'warnings':list(event.warnings),'differences':[],'segments':[],
            'reference_match':'user_confirmed' if confirmed else 'unconfirmed',
            'reference_assignment_confirmed':bool(confirmed),
            'actual_sample_identity_status':'NOT_ESTABLISHED',
            'coordinate_convention':'0-based half-open; query coordinates use original FASTA; query_oriented coordinates use aligned strand',
            'native_snapgene_alignment':False}
    if not event.consensus:return result
    alignment=align_consensus(reference,event.consensus);result.update({k:alignment[k] for k in ('orientation','differences','segments','matched_bases','query_length','reference_rotation0')})
    cov=alignment['coverage'];mapping=alignment['query_mapping'];bm=event.bam_metrics
    quality=bool(bm.get('depth')) and bm.get('reference_length')==len(event.consensus) and bm.get('forward',0)>0 and bm.get('reverse',0)>0
    depth=bm.get('depth',[])
    def supported(pos):
        qp=mapping.get(pos)
        if not quality or qp is None:return cov.get(pos) in 'ACGT'
        strands_ok=all(values[qp]>=1 for values in (bm.get('depth_forward',depth),bm.get('depth_reverse',depth)))
        return cov.get(pos) in 'ACGT' and depth[qp]>=minimum_depth and strands_ok
    def complete(reg):return all(pos in cov and supported(pos) for pos in range(reg.start,reg.end))
    def deletion_support(d,anchor=8):
        left=[(d['start0']-anchor+i)%len(reference.sequence) for i in range(anchor)]
        right=[(d['start0']+d['length']+i)%len(reference.sequence) for i in range(anchor)]
        positions=left+right
        proof={'anchor_bp':anchor,'unique_spanning_reads':0,'forward_spanning_reads':0,'reverse_spanning_reads':0,'supported':False}
        if not quality or not all(p in mapping and p in cov and supported(p) and cov[p]==reference.sequence[p] for p in positions):
            proof['reason']='INSUFFICIENT_MAPPED_BIDIRECTIONAL_ANCHORS';return proof
        qp=[mapping[p] for p in positions];lo,hi=min(qp),max(qp)+1
        if hi-lo!=2*anchor or len(set(qp))!=2*anchor:
            proof['reason']='ANCHORS_DO_NOT_BRIDGE_ONE_QUERY_BOUNDARY';return proof
        names=[set(),set()]
        for block in bm.get('read_blocks',[]):
            if block['start0']<=lo and block['end0']>=hi:
                offset=lo-block['start0']
                if block['sequence'][offset:offset+hi-lo]==event.consensus[lo:hi]:
                    names[bool(block['reverse'])].add(block['read_name'])
        proof.update({'query_anchor_start0':lo,'query_anchor_end0':hi,'unique_spanning_reads':len(names[0]|names[1]),
                      'forward_spanning_reads':len(names[0]),'reverse_spanning_reads':len(names[1])})
        proof['supported']=len(names[0]|names[1])>=minimum_variant_reads and bool(names[0]) and bool(names[1])
        proof['reason']='SUPPORTED_DELETION_JUNCTION' if proof['supported'] else 'INSUFFICIENT_ACTUAL_SPANNING_READS'
        return proof
    for d in alignment['differences']:
        if d['type']=='DEL':d['support_evidence']=deletion_support(d)
    for reg in reference.regions:
        diffs=[d for d in alignment['differences'] if reg.name in d['required_regions']]
        coverage=sum(pos in cov for pos in range(reg.start,reg.end))/(reg.end-reg.start)
        confirmed_differences=[]
        for d in diffs:
            if d['type']=='DEL':
                if d['support_evidence']['supported']:confirmed_differences.append(d)
            elif d['type']=='INS':
                pos=d['start0'];previous=(pos-1)%len(reference.sequence)
                if quality and pos in mapping and previous in mapping and supported(pos) and supported(previous) and 0<d['query_start0']<len(event.consensus)-d['length']:
                    confirmed_differences.append(d)
            elif d['alt'] in 'ACGT':
                pos=d['start0']
                if quality and pos in mapping and pos in cov and supported(pos):confirmed_differences.append(d)
        status='FAIL' if confirmed_differences else 'INCOMPLETE' if diffs else 'PASS' if complete(reg) and quality else 'INCOMPLETE'
        result['regions'].append({'region':reg.name,'kind':reg.kind,'start0':reg.start,'end0':reg.end,'required':reg.required,
                                  'coverage_fraction':coverage,'status':status})
    failures=any(x['status']=='FAIL' and x['required'] for x in result['regions'])
    incomplete=any(x['status']=='INCOMPLETE' and x['required'] for x in result['regions'])
    result['target_status']='FAIL' if failures else 'INCOMPLETE' if incomplete else 'PASS'
    # Empty vendor tables never validate the design. Preserve all heterogeneity evidence.
    result['vendor_variant_evidence']=event.variants
    mixtures=[]
    for d in bm.get('large_deletions',[]):
        if d['unique_reads']>=minimum_variant_reads and d['fraction']>=.05:
            mixtures.append({'source':'BAM_CIGAR','type':'DEL',**d})
    inverse_mapping={qp:rp for rp,qp in mapping.items()}
    for row in event.variants:
        try:
            candidate=float(row.get('Freq',0))>=.05 and int(row.get('Alt_dp',0))>=minimum_variant_reads
            query_pos=int(row.get('Pos',0))-1
        except (ValueError,TypeError):candidate=False;query_pos=-1
        if row.get('table')=='filtered' or candidate:
            design_pos=inverse_mapping.get(query_pos)
            mixtures.append({'source':'vendor_retained_variant' if row.get('table')=='filtered' else 'vendor_unfiltered_candidate_requires_review',
                **row,'design_position0':design_pos,'required_regions':[reg.name for reg in reference.regions if design_pos is not None and reg.start<=design_pos<reg.end]})
    result['mixture_evidence']=mixtures
    backbone=[d for d in alignment['differences'] if not d['required_regions']]
    if mixtures and result['target_status']!='FAIL':
        result['target_status']='REVIEW'
        if incomplete:result['warnings'].append('MIXTURE_AND_REQUIRED_COVERAGE_INCOMPLETE')
    elif (backbone or 'AB1_CONSENSUS_CONFLICT' in event.warnings) and result['target_status']=='PASS':result['target_status']='REVIEW'
    result['expected_design_status']=result['target_status']
    result['expected_design_mismatch_confirmed']=bool(confirmed and failures)
    if not confirmed:
        result['computed_target_status']=result['target_status'];result['target_status']='REVIEW';result['warnings'].append('REFERENCE_ASSIGNMENT_NOT_CONFIRMED')
    if 'SAMPLE_NAME_COLLISION' in event.warnings:
        result['warnings'].append('CLONE_RELATION_UNKNOWN_UNLESS_EXPLICITLY_RECORDED')
        if not event.colony_id and result['target_status']=='PASS':result['target_status']='REVIEW'
    if len(cov)==len(reference.sequence) and quality and all(supported(pos) for pos in range(len(reference.sequence))):
        result['whole_plasmid_status']='FAIL' if alignment['differences'] else 'REVIEW' if mixtures else 'PASS'
        if not confirmed:result['whole_plasmid_status']='REVIEW'
    if not quality:result['warnings'].append('INDEPENDENT_READ_SUPPORT_INCOMPLETE')
    event.construct_identity=reference.construct;event.reference_sha256=reference.sequence_sha256
    return result

def reference_candidates(event,references):
    """Names generate candidates only. No automatic identity assignment."""
    label=event.sample_name.rsplit('-',1)[0]
    found=[]
    for ref in references:
        exact=label in (ref.construct,ref.operational_label,*ref.aliases)
        km={event.consensus[i:i+21] for i in range(19,max(19,len(event.consensus)-20),21)}
        hits=sum(k in ref.sequence+ref.sequence[:len(event.consensus)] for k in km)
        if exact or hits:found.append({'construct':ref.construct,'sha256':ref.sequence_sha256,'name_evidence':exact,'shared_21mers':hits})
    return sorted(found,key=lambda x:(x['name_evidence'],x['shared_21mers']),reverse=True)


def unassigned_result(event):
    return {'event_id':event.event_id,'company_id':event.company_id,'sample_name':event.sample_name,
        'batch':event.processing_batch,'clone_id':event.clone_id,'colony_id':event.colony_id,
        'colony_suffix':event.colony_suffix,'construct':'UNASSIGNED','reference_sha256':'',
        'reference_snapshot_sha256':'','policy_sha256':'','consensus_sha256':event.consensus_sha256,
        'target_status':'REVIEW' if event.consensus else 'INCOMPLETE',
        'whole_plasmid_status':'INCOMPLETE','regions':[],'differences':[],'segments':[],
        'warnings':list(event.warnings)+['REFERENCE_NOT_ASSIGNED'],'native_snapgene_alignment':False}
