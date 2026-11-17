"""Separate user-scoped verification ledger. Never opens the production primer database."""
from dataclasses import asdict
import json,sqlite3,hashlib
from pathlib import Path
from verification_import import Event,merge_events

def load_events(path):
    p=Path(path)
    if not p.exists():return []
    with sqlite3.connect('file:'+str(p.resolve())+'?mode=ro',uri=True) as c:
        return [Event(**json.loads(row[0])) for row in c.execute('select payload from sequencing_events order by rowid')]

def save_verifications(path,events,references,results):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True,mode=0o700)
    if p.is_symlink() or p.name in ('app.db','primer.db'):raise ValueError('Production database path forbidden')
    existing=load_events(p) if p.exists() else []
    merged=merge_events(existing,events)
    with sqlite3.connect(p) as c:
        c.execute('create table if not exists sequencing_events(event_id text primary key, payload text not null)')
        c.execute('create table if not exists frozen_references(construct text, sha256 text, payload text, primary key(construct,sha256))')
        c.execute('create table if not exists verifications(event_id text, reference_sha256 text, payload text, primary key(event_id,reference_sha256))')
        c.execute('create table if not exists verification_runs(run_sha256 text primary key, event_id text, reference_sha256 text, policy_sha256 text, payload text)')
        for e in merged:c.execute('insert or replace into sequencing_events values(?,?)',(e.event_id,json.dumps(asdict(e))))
        c.execute('create table if not exists reference_snapshots(snapshot_sha256 text primary key,payload text)')
        for ref in references:
            payload=json.dumps(ref.snapshot(),sort_keys=True)
            c.execute('insert or ignore into reference_snapshots values(?,?)',(hashlib.sha256(payload.encode()).hexdigest(),payload))
            c.execute('insert or ignore into frozen_references values(?,?,?)',(ref.construct,ref.sequence_sha256,payload))
        for v in results:
            payload=json.dumps(v,sort_keys=True)
            c.execute('insert or ignore into verification_runs values(?,?,?,?,?)',(hashlib.sha256(payload.encode()).hexdigest(),v['event_id'],v['reference_sha256'],v.get('policy_sha256',''),payload))
            c.execute('insert or replace into verifications values(?,?,?)',(v['event_id'],v['reference_sha256'],payload))
    p.chmod(0o600)


def load_ledger(path):
    from verification_reference import Reference,Region
    p=Path(path)
    if not p.exists():return [],[],[]
    events=load_events(p);references=[];results=[]
    with sqlite3.connect('file:'+str(p.resolve())+'?mode=ro',uri=True) as c:
        tables={x[0] for x in c.execute("select name from sqlite_master where type='table'")}
        if 'reference_snapshots' in tables:
            for snapshot_hash,payload in c.execute('select snapshot_sha256,payload from reference_snapshots'):
                if hashlib.sha256(payload.encode()).hexdigest()!=snapshot_hash:raise ValueError('Stored snapshot checksum mismatch')
                item=json.loads(payload);item['regions']=tuple(Region(**reg) for reg in item['regions']);item['aliases']=tuple(item['aliases'])
                ref=Reference(**item)
                if hashlib.sha256(ref.sequence.encode()).hexdigest()!=ref.sequence_sha256:raise ValueError('Stored reference checksum mismatch')
                references.append(ref)
        if 'verification_runs' in tables:
            results=[json.loads(row[0]) for row in c.execute('select payload from verification_runs order by rowid')]
    return events,references,results
