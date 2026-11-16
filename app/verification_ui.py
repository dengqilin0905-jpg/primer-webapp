"""Sequencing verification UI: explicit assignment, independent private ledger and exports."""
from pathlib import Path
import io,json,os,hashlib
from dataclasses import asdict
from Bio import SeqIO
import streamlit as st
from verification_import import parse_zip,merge_events,Event
from verification_reference import freeze_reference,require_features
from verification_catalog import load_catalog
from verification_align import reference_candidates,verify_event,unassigned_result
from verification_export import project_exports
from verification_store import load_ledger,save_verifications
from project_contract import PROJECTS_ROOT

def sequencing_verification_page(user,cfg):
    st.subheader('Sequencing Verification')
    st.caption('独立设计参考验证；公司自身 consensus 比对及生成式 AB1 不代表构建通过。Target 与全质粒状态独立。FAIL 说明与已确认的预期设计不符，不代表实际样品身份已查明。')
    uid=str(int(user['id']));prefix='verification_'+uid
    events=st.session_state.setdefault(prefix+'_events',[])
    refs=st.session_state.setdefault(prefix+'_refs',[])
    ledger=Path(os.environ.get('VERIFICATION_STATE_DIR',str(Path.home()/'.primer-verification')))/uid/'verification.sqlite3'
    if st.button('载入独立私有验证历史',key=prefix+'_load'):
        try:
            saved_events,saved_refs,saved_results=load_ledger(ledger)
            events=merge_events(events,saved_events);st.session_state[prefix+'_events']=events
            known={(x.construct,x.sequence_sha256) for x in refs}
            for ref in saved_refs:
                if (ref.construct,ref.sequence_sha256) not in known:
                    refs.append(ref);known.add((ref.construct,ref.sequence_sha256))
            st.session_state[prefix+'_history_results']=saved_results
        except Exception as exc:st.error(f'历史未载入：{exc}')
    if st.button('读取服务器正式项目参考（只读）',key=prefix+'_catalog'):
        refs,issues=load_catalog(PROJECTS_ROOT,Path(__file__).resolve().parents[1]);st.session_state[prefix+'_refs']=refs
        st.info(f'已冻结 {len(refs)} 条可验证参考；{len(issues)} 条需要补充 provenance。')
        if issues:st.text('\n'.join(issues))
    with st.expander('上传参考 GBK 及明确组装 provenance'):
        gb=st.file_uploader('正式参考',type=['gb','gbk'],key=prefix+'_gb')
        pv=st.file_uploader('组装 JSON：backbone_sequence, replace_start/end, insert_sequence, source, operational_label',type=['json'],key=prefix+'_provenance')
        if st.button('冻结上传参考',key=prefix+'_freeze') and gb and pv:
            try:
                record=SeqIO.read(io.StringIO(gb.getvalue().decode()),'genbank');data=json.loads(pv.getvalue())
                ref=freeze_reference(record,**data)
                if not any(x.construct==ref.construct and x.sequence_sha256==ref.sequence_sha256 for x in refs):refs.append(ref)
            except Exception as exc:st.error(f'参考未冻结：{exc}')
    uploads=st.file_uploader('康为结果 ZIP（可跨批次多选）',type=['zip'],accept_multiple_files=True,key=prefix+'_zips')
    if st.button('导入结果包',key=prefix+'_import'):
        for upload in uploads:
            try:events=merge_events(events,parse_zip(upload.getvalue()))
            except Exception as exc:st.error(f'{upload.name} 未导入：{exc}')
        st.session_state[prefix+'_events']=events
    history=st.file_uploader('恢复此前下载的 events.json',type=['json'],key=prefix+'_json')
    if st.button('恢复 JSON 历史',key=prefix+'_restore') and history:
        try:
            rows=json.loads(history.getvalue())
            if len(rows)>200:raise ValueError('History row limit')
            events=merge_events(events,[Event(**row) for row in rows]);st.session_state[prefix+'_events']=events
        except Exception as exc:st.error(f'历史未恢复：{exc}')
    st.dataframe([{'公司编号':e.company_id,'批次':e.processing_batch,'样品':e.sample_name,'colony后缀':e.colony_suffix,'物理colony':e.colony_id or 'unknown','事件':e.event_id[:12],'警告':';'.join(e.warnings)} for e in events],use_container_width=True)
    historical=st.session_state.get(prefix+'_history_results',[])
    if historical:
        batch=st.selectbox('历史批次',['全部']+sorted({x.get('batch','') for x in historical}),key=prefix+'_filter_batch')
        construct=st.selectbox('历史构建',['全部']+sorted({x.get('construct','') for x in historical}),key=prefix+'_filter_construct')
        colony=st.text_input('历史 colony 筛选（ID或后缀）',key=prefix+'_filter_colony').strip()
        historical=[x for x in historical if (batch=='全部' or x.get('batch')==batch) and (construct=='全部' or x.get('construct')==construct) and (not colony or colony in str(x.get('colony_id','')) or colony==str(x.get('colony_suffix','')))]
        with st.expander('已保存的跨批次验证记录'):
            st.dataframe([{k:x.get(k,'') for k in ('company_id','sample_name','batch','clone_id','colony_id','construct','target_status','whole_plasmid_status','policy_sha256')} for x in historical],use_container_width=True)
            for x in historical:
                with st.container():
                    st.caption(x['company_id']+' / '+x['reference_sha256'][:12]+' / '+x['policy_sha256'][:12])
                    st.dataframe(x['regions']);st.dataframe(x['differences'])
    assignments={}
    for e in events:
        with st.expander(e.company_id+' / '+e.sample_name+' / '+e.processing_batch):
            st.dataframe(reference_candidates(e,refs),use_container_width=True)
            options=['不匹配']+[ref.construct+' / '+ref.sequence_sha256[:12] for ref in refs]
            choice=st.selectbox('正式参考（候选仅供核查，不自动认定）',options,key=prefix+e.event_id+'_ref')
            e.clone_id=st.text_input('克隆轮次 / clone ID（不知道则留空）',value=e.clone_id,key=prefix+e.event_id+'_clone')
            e.colony_id=st.text_input('物理 colony ID（后缀不等于唯一身份）',value=e.colony_id,key=prefix+e.event_id+'_colony')
            confirmed=st.checkbox('已核对送检预期参考对应关系（不代表已查明样品实际身份）',key=prefix+e.event_id+'_confirmed')
            if choice!='不匹配':
                ref=refs[options.index(choice)-1]
                record=SeqIO.read(io.StringIO(ref.genbank),'genbank')
                labels=sorted({';'.join(f.qualifiers.get('label',[])) for f in record.features if f.qualifiers.get('label')})
                required=st.multiselect('Additional required functional elements',labels,key=prefix+e.event_id+'_elements')
                assignments[e.event_id]=(require_features(ref,required),confirmed)
    token=hashlib.sha256(json.dumps({'events':[e.event_id for e in events],'assignments':[(e.event_id,e.clone_id,e.colony_id,assignments[e.event_id][0].sequence_sha256,assignments[e.event_id][0].provenance_json,assignments[e.event_id][1]) for e in events if e.event_id in assignments]},sort_keys=True).encode()).hexdigest()
    if st.button('运行独立序列验证',key=prefix+'_run'):
        results=[]
        for e in events:
            if e.event_id not in assignments:
                results.append(unassigned_result(e));continue
            ref,confirmed=assignments[e.event_id]
            try:results.append(verify_event(e,ref,confirmed=confirmed))
            except Exception as exc:
                st.error(f'{e.company_id} verification incomplete: {exc}')
                failed=unassigned_result(e);failed['target_status']='INCOMPLETE';failed['construct']=ref.construct
                failed['reference_sha256']=ref.sequence_sha256
                failed['reference_snapshot_sha256']=hashlib.sha256(json.dumps(ref.snapshot(),sort_keys=True).encode()).hexdigest()
                failed['warnings'].append('ALIGNMENT_NOT_COMPLETED');results.append(failed)
        st.session_state[prefix+'_results']=results
        st.session_state[prefix+'_result_token']=token
    results=st.session_state.get(prefix+'_results',[])
    if st.session_state.get(prefix+'_result_token')!=token:
        results=[]
    if results:
        st.dataframe([{k:x[k] for k in ('company_id','sample_name','batch','colony_id','construct','target_status','whole_plasmid_status')} for x in results],use_container_width=True)
        for x in results:
            with st.expander('QC / '+x['company_id']):
                st.dataframe(x['regions']);st.dataframe(x['differences']);st.json({'warnings':x['warnings'],'mixtures':x.get('mixture_evidence',[])})
        export_refs=[assignments[x['event_id']][0] for x in results if x['event_id'] in assignments]
        outputs=project_exports(results,export_refs,events)
        for kind,mime in [('xlsx','application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'),('html','text/html'),('zip','application/zip')]:
            st.download_button('下载验证结果 '+kind,outputs[kind],file_name='sequencing_verification.'+kind,mime=mime,key=prefix+'_download_'+kind)
    if events:
        st.download_button('下载跨批次事件历史 JSON',json.dumps([asdict(e) for e in events],ensure_ascii=False).encode(),file_name='events.json',key=prefix+'_history')
        if st.button('保存到独立私有验证库',key=prefix+'_save'):
            try:save_verifications(ledger,events,[a[0] for a in assignments.values()],results);st.success('已保存独立验证记录；未修改 primer 生产数据库或实物库存。')
            except Exception as exc:st.error(f'未保存：{exc}')
