import hashlib, json
"""UI adapter for the existing exporter; no primer generation or historical rewrites."""
from Bio.Seq import Seq
from genbank_export import build_genbank_record
from naming import build_operational_label
from sequencing_export import export_sequencing_workbook,resolve_sequencing_template,antibiotic_options
from openpyxl import load_workbook
import io
from identity_audit import load_identity_tsv,identity_key,RICE_IDENTITY_REGISTRY
from project_contract import canonical_species_abbreviation

def items_from_design(rows,metadata):
    if len(rows)!=len(metadata):raise ValueError('Metadata row count mismatch')
    result=[]
    for row,meta in zip(rows,metadata):
        hf=str(row.get('homo_f') or '').upper()
        if len(hf)!=18 or set(hf)-set('ACGT'):
            raise ValueError('需要已确认的 18 bp 同源臂；历史 21 bp 不自动截短，请按原图谱确认有效 profile。')
        record=build_genbank_record(row)
        if str(record.seq).count(hf)!=1:raise ValueError('Reference arm must locate uniquely in actual map')
        symbol=str(meta.get('symbol') or '').strip();species=str(meta.get('species') or '').strip()
        identity=str(meta.get('canonical_identity') or '').strip()
        if identity and canonical_species_abbreviation(species)=='Os':
            gene=str(meta.get('gene_id') or '').strip()
            evidence=load_identity_tsv(RICE_IDENTITY_REGISTRY,source_kind='identity_registry')
            matches={identity_key(x.canonical_identity) for x in evidence if x.gene_id==gene}
            if matches!={identity_key(identity)}:
                raise ValueError('正式Rice identity与Gene ID不一致，或registry证据缺失。')
        label=build_operational_label(symbol,species=species or None,
              identifier=str(meta.get('canonical_identity') or '').strip() or None,
              backbone=str(row.get('backbone_label') or row.get('backbone_name') or ''))
        resistance=str(meta.get('resistance') or '').strip()
        if not resistance:raise ValueError('请根据真实图谱确认抗性，不能从载体名称猜测。')
        result.append({'operational_label':label,'full_construct_name':record.id,
            'reference_primer_name':label+'_F','reference_sequence':hf,'plasmid_length':len(record),
            'sequencing_length':min(len(row['insert_seq'])+18,1499),'resistance':resistance,
            'source_map':record.id+'.gbk','remarks':''})
    return result

def sequencing_download(rows,key):
    import streamlit as st
    with st.expander('康为四代送测单',expanded=False):
        st.caption('复用已完成设计，不重新设计引物。确认正式 species、identity/symbol 及图谱抗性；默认每构建 3 个 colony。软件内置为通用示例模板，正式送测请配置当前公司原模板。')
        defaults=[{'symbol':str(row.get('gene_name') or ''),'species':str(row.get('species') or ''),
                   'canonical_identity':str(row.get('canonical_identity') or ''),'gene_id':str(row.get('locus_tag') or ''),'resistance':''} for row in rows]
        options=antibiotic_options(load_workbook(io.BytesIO(resolve_sequencing_template().template_bytes)))
        metadata=st.data_editor(defaults,key=key+'_metadata',num_rows='fixed',column_config={'resistance':st.column_config.SelectboxColumn('抗性（模板正式值）',options=list(options))})
        confirmed=st.checkbox('已按正式 identity 记录及图谱确认以上字段',key=key+'_confirmed')
        token=hashlib.sha256(json.dumps([rows,metadata],sort_keys=True,default=str).encode()).hexdigest()
        if st.button('生成康为送测 Excel',key=key+'_generate',disabled=not confirmed):
            try:
                data,report=export_sequencing_workbook(items_from_design(rows,metadata))
                st.session_state[key+'_artifact']=(data,report,token)
            except Exception as exc:st.error(f'送测单未生成：{exc}')
        if key+'_artifact' in st.session_state and st.session_state[key+'_artifact'][2]==token and confirmed:
            data,report,_=st.session_state[key+'_artifact']
            st.download_button('下载康为送测 Excel',data,file_name='sequencing_order.xlsx',key=key+'_download')
            st.caption(f"样本数：{report['sample_count']}；模板：{report['template_source']}")
