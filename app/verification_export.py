"""Self-contained review exports; never claim native SnapGene alignment output."""
from dataclasses import asdict
from html import escape
from hashlib import sha256
import io,json,csv,zipfile,re
from openpyxl import Workbook
from verification_html import render_alignment

def table_tsv(rows):
    if not rows:return ''
    fields=list(dict.fromkeys(k for row in rows for k in row));s=io.StringIO();w=csv.DictWriter(s,fieldnames=fields,delimiter='\t');w.writeheader()
    for row in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()})
    return s.getvalue()

def project_exports(results,references,events):
    refs={sha256(json.dumps(x.snapshot(),sort_keys=True).encode()).hexdigest():x for x in references};evs={e.event_id:e for e in events}
    wb=Workbook();wb.remove(wb.active)
    sheets={'Samples':[{k:v for k,v in x.items() if k not in ('regions','differences','segments','mixture_evidence')} for x in results],
            'Regions':[{'event_id':x['event_id'],**row} for x in results for row in x['regions']],
            'Differences':[{'event_id':x['event_id'],**row} for x in results for row in x['differences']],
            'Evidence':[{'event_id':x['event_id'],'mixtures':x.get('mixture_evidence',[]),'vendor_candidates':x.get('vendor_variant_evidence',[]),'warnings':x['warnings']} for x in results]}
    for name,rows in sheets.items():
        ws=wb.create_sheet(name);fields=list(dict.fromkeys(k for row in rows for k in row));ws.append(fields)
        for row in rows:
            ws.append([json.dumps(row.get(k,''),ensure_ascii=False) if isinstance(row.get(k),(dict,list)) else row.get(k,'') for k in fields])
            for cell in ws[ws.max_row]:
                if isinstance(cell.value,str):cell.data_type='s'
        ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions
    bio=io.BytesIO();wb.save(bio)
    html=['<!doctype html><html lang="zh"><meta charset="utf-8"><title>Sequencing Verification</title><style>body{font:15px sans-serif;margin:32px;color:#172b3a}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ccd8df;padding:8px;text-align:left}pre{white-space:pre-wrap;overflow-wrap:anywhere}.base-alignment{margin:20px 0}.alignment-block{overflow-x:auto;margin:10px 0;border:1px solid #dce4ea;padding:8px}.alignment-row{white-space:pre!important;overflow-wrap:normal!important;font:13px monospace;margin:3px 0}.line-label{display:inline-block;width:260px;color:#345}.base{display:inline-block;width:1ch}.mismatch{background:#ffd2d2;color:#8a1010}.insertion{background:#c7f1d8;color:#075327}.deletion{background:#ffe4b5;color:#794500}.feature-label{font:12px sans-serif;color:#52646e}.difference-context{margin:12px 0}.homopolymer{font-weight:bold}</style><h1>Sequencing Verification</h1><p>Target and whole-plasmid statuses are independent. No native SnapGene alignment file is produced. Exported coordinates are 0-based half-open; base views use 1-based original-sequence coordinates.</p><table><tr><th>Sample</th><th>Company / batch</th><th>Construct</th><th>Target</th><th>Whole plasmid</th></tr>']
    for x in results:html.append('<tr>'+''.join('<td>'+escape(str(v))+'</td>' for v in (x['sample_name'],x['company_id']+' / '+x['batch'],x['construct'],x['target_status'],x['whole_plasmid_status']))+'</tr>')
    html.append('</table>')
    for x in results:
        ref=refs.get(x['reference_snapshot_sha256']);e=evs[x['event_id']]
        html.append('<article><h2>'+escape(x['company_id']+' '+x['sample_name'])+'</h2>')
        if ref is not None:html.append(render_alignment(x,ref,e))
        html.append('<details><summary>QC evidence / JSON</summary><pre>'+escape(json.dumps(x,ensure_ascii=False,indent=2))+'</pre></details></article>')
    html.append('</html>');html=''.join(html).encode()
    out=io.BytesIO()
    with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
        z.writestr('verification.xlsx',bio.getvalue());z.writestr('verification.html',html)
        z.writestr('events.json',json.dumps([asdict(e) for e in events],ensure_ascii=False,indent=2))
        for x in results:
            ref=refs.get(x['reference_snapshot_sha256']);e=evs[x['event_id']]
            folder=re.sub(r'[^A-Za-z0-9_.-]','_',e.company_id+'_'+e.sample_name)+'_'+e.event_id[:12]
            if ref is not None:
                z.writestr(folder+'/reference.gbk',ref.genbank);z.writestr(folder+'/reference_snapshot.json',json.dumps(ref.snapshot(),ensure_ascii=False,indent=2))
            z.writestr(folder+'/consensus.fasta','>'+e.company_id+' '+e.sample_name+'\n'+e.consensus+'\n')
            if ref is not None:z.writestr(folder+'/reference.fasta','>'+ref.construct+'\n'+ref.sequence+'\n')
            z.writestr(folder+'/alignment_coordinates.tsv',table_tsv(x['segments']))
            z.writestr(folder+'/differences.tsv',table_tsv(x['differences']))
            z.writestr(folder+'/qc.json',json.dumps(x,ensure_ascii=False,indent=2))
    return {'xlsx':bio.getvalue(),'html':html,'zip':out.getvalue()}
