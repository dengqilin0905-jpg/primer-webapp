"""Actual aligned bases reconstructed from the exported alignment segments, not JSON summaries."""
from html import escape
from hashlib import sha256
import io
from Bio import SeqIO
from Bio.Seq import Seq


def alignment_columns(result,reference,event):
    if reference.sequence_sha256!=result['reference_sha256']:
        raise ValueError('Alignment reference version mismatch')
    if event.consensus_sha256 and sha256(event.consensus.encode()).hexdigest()!=event.consensus_sha256:
        raise ValueError('Consensus checksum mismatch')
    raw=event.consensus;direction=result.get('orientation','+')
    oriented=raw if direction=='+' else str(Seq(raw).reverse_complement())
    rec=SeqIO.read(io.StringIO(reference.genbank),'genbank');cols=[];n=len(reference.sequence)
    for segment in result['segments']:
        a=segment['reference_start0'];dx=segment['reference_unwrapped_end0']-segment['reference_unwrapped_start0']
        # Legacy v1 segments used oriented coordinates in query_start0/end0.
        b=segment.get('query_oriented_start0',segment['query_start0'])
        e=segment.get('query_oriented_end0',segment['query_end0']);dy=e-b
        if dx and dy and dx!=dy:raise ValueError('Unsupported non-diagonal alignment segment')
        for j in range(max(dx,dy)):
            rp=(a+j)%n if dx else None;oq=b+j if dy else None
            qp=oq if direction=='+' else len(raw)-1-oq if oq is not None else None
            rb=reference.sequence[rp] if rp is not None else '-';qb=oriented[oq] if oq is not None else '-'
            anchor=rp if rp is not None else a
            required=[x.name for x in reference.regions if x.required and x.start<=anchor<x.end]
            features=[]
            for f in rec.features:
                if f.type!='source' and len(f)<n and anchor in f.location:
                    features.extend(f.qualifiers.get('label',[f.type]))
            scope='required_target' if required else 'retained_backbone'
            kind='insertion' if rp is None else 'deletion' if oq is None else 'mismatch' if rb!=qb else 'match'
            cols.append({'ref':rb,'query':qb,'reference_pos0':rp,'query_pos0':qp,'query_oriented_pos0':oq,
                         'anchor_reference_pos0':anchor,'kind':kind,'scope':scope,'features':features,'required_regions':required})
    return cols


def _range_label(values):
    v=[x+1 for x in values if x is not None]
    return f'{v[0]}→{v[-1]}' if v else 'gap'


def _block(columns,direction):
    labels=sorted({x for c in columns for x in c['features']})
    scopes=sorted({c['scope'] for c in columns})
    rows=[]
    for field,prefix,pos in [('ref','Reference +','reference_pos0'),('query','Query '+direction+' original','query_pos0')]:
        cells=[]
        for c in columns:
            title='; '.join([c['kind'],c['scope'],*c['features']])
            attrs=f'data-ref0="{c["reference_pos0"] if c["reference_pos0"] is not None else ""}" data-query0="{c["query_pos0"] if c["query_pos0"] is not None else ""}"'
            cells.append(f'<span class="base {c["kind"]}" {attrs} title="{escape(title,quote=True)}">{c[field]}</span>')
        rows.append('<pre class="alignment-row"><span class="line-label">'+escape(prefix+' '+_range_label([c[pos] for c in columns]))+'</span>'+''.join(cells)+'</pre>')
    return '<div class="alignment-block"><div class="feature-label">'+escape(' / '.join(scopes)+'; features: '+(', '.join(labels) or 'unannotated backbone'))+'</div>'+''.join(rows)+'</div>'


def render_alignment(result,reference,event):
    if not event.consensus:return '<p>No sequence returned; base alignment unavailable (INCOMPLETE).</p>'
    if not result['segments']:return '<p>Alignment not completed; no base-level validation is claimed.</p>'
    cols=alignment_columns(result,reference,event);direction=result['orientation']
    parts=['<section class="base-alignment"><h3>实际碱基比对 / Actual aligned bases</h3>',
           '<p>显示坐标为 1-based；下载 TSV 为 0-based half-open。Reference 按原质粒正向坐标；Query 使用原始 FASTA 坐标，反向比对时递减；定向后坐标另列于 TSV。“-”表示 gap。</p>',
           '<p class="legend"><span class="mismatch">mismatch</span> <span class="insertion">insertion</span> <span class="deletion">deletion / gap</span></p>',
           '<p>长缺失、低相似度和重复序列的 gap 位置可能不唯一；此视图不是独立读段确认，也不说明实际样品身份已查明。确认缺失还需两端锚点及实际跨接读段证据。</p>']
    # Important contexts are shown first, followed by the complete returned alignment.
    parts.append('<h4>差异上下文 / Difference contexts</h4>')
    if not result['differences']:parts.append('<p class="no-differences">No differences; all returned aligned bases agree.</p>')
    for i,d in enumerate(result['differences']):
        indices=[j for j,c in enumerate(cols) if (d['type']=='INS' and c['kind']=='insertion' and c['anchor_reference_pos0']==d['start0']) or (d['type']!='INS' and c['reference_pos0'] is not None and (c['reference_pos0']-d['start0'])%len(reference.sequence)<(d['length'] if d['type']=='DEL' else 1))]
        if not indices:continue
        lo,hi=max(0,min(indices)-20),min(len(cols),max(indices)+21)
        label=f'{d["type"]} @ reference {d["start0"]+1}; {d["classification"]}; features: '+(', '.join(d['features']) or 'unannotated backbone')
        parts.append('<details class="difference-context" open><summary>'+escape(label)+'</summary>')
        if d['type']=='INS' and d['alt'] and len(set(d['alt']))==1:
            anchor=d['start0'];s=reference.sequence;base=d['alt'][0];left=anchor;right=anchor
            while left>0 and s[left-1]==base:left-=1
            while right<len(s) and s[right]==base:right+=1
            if right-left>=4:parts.append('<p class="homopolymer">'+escape(f'{base}{right-left} → {base}{right-left+len(d["alt"])}; homopolymer placement is nonunique; '+d['classification'])+'</p>')
        if d['type']=='DEL':parts.append('<p>'+escape('Deletion support: '+str(d.get('support_evidence','not independently evaluated')) )+'</p>')
        for j in range(lo,hi,70):parts.append(_block(cols[j:min(j+70,hi)],direction))
        parts.append('</details>')
    parts.append('<details class="complete-alignment"><summary>完整返回片段比对 / Complete returned alignment</summary>')
    for j in range(0,len(cols),70):parts.append(_block(cols[j:j+70],direction))
    parts.append('</details></section>')
    return ''.join(parts)
