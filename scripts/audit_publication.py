"""Audit software snapshot and new Git history without printing matched secrets."""
import argparse,subprocess,re,sys,zipfile,io
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ALLOWED_TOP={'PUBLICATION_FILES.txt','.env.server.example','.gitignore','AGENT.md','CHANGELOG.md','DEPLOYMENT.md','LICENSE','README.md','README.txt','WORK_LOG.md','requirements.txt','requirements-lock.txt','pyproject.toml','__init__.py','launcher.py','run.bat','run_debug.bat'}
ALLOWED_DIRS={'.github','app','config','deploy','docs','examples','resources','scripts','test'}
PROHIBITED={'.dna','.bam','.bai','.ab1','.fastq','.fq','.db','.sqlite','.sqlite3','.zip'}
PATTERNS={
 'personal_home':rb'/home/[A-Za-z0-9_.-]+/',
 'personal_windows':rb'[A-Za-z]:\\Users\\',
 'company_record':rb'S2[0-9]{10,}',
 'private_key':rb'-----BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY-----',
 'github_token':rb'(?:gh[pousr]_[A-Za-z0-9]{25,}|github_pat_[A-Za-z0-9_]{30,})',
 'cloud_key':rb'AKIA[0-9A-Z]{16}',
}

def inspect(path,data):
    p=Path(path);issues=[]
    if p.suffix.lower() in PROHIBITED:issues.append('private-data file type')
    if p.suffix.lower() in {'.fa','.fasta','.gb','.gbk'} and not (p.parts[0]=='examples' and p.name.startswith('synthetic_')):issues.append('non-synthetic sequence asset')
    if p.suffix.lower()=='.xlsx':
        if p.parent.as_posix()!='resources/templates' or not p.name.endswith('_generic.xlsx'):issues.append('unreviewed workbook')
        else:
            with zipfile.ZipFile(io.BytesIO(data)) as z:
                data=b'\n'.join(z.read(n) for n in z.namelist() if n.endswith('.xml'))
    issues.extend(name for name,pattern in PATTERNS.items() if re.search(pattern,data))
    return issues

def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--base',default='origin/main');args=parser.parse_args();failures=[]
    paths=set(git('ls-files','--cached','--others','--exclude-standard').decode().splitlines())
    for path in sorted(paths):
        p=ROOT/path
        if not p.is_file():continue
        if p.is_symlink():failures.append((path,['symlink not in software allowlist']));continue
        if Path(path).parts[0] not in ALLOWED_DIRS and path not in ALLOWED_TOP:failures.append((path,['not in software allowlist']));continue
        found=inspect(path,p.read_bytes())
        if found:failures.append((path,found))
    new=git('rev-list','--objects',args.base+'..HEAD').decode().splitlines();blobs=0
    for line in new:
        oid,_,path=line.partition(' ')
        if git('cat-file','-t',oid).strip()!=b'blob':continue
        blobs+=1;found=inspect(path,git('cat-file','blob',oid))
        if found:failures.append(('history:'+oid[:12]+':'+path,found))
    for path,issues in failures:print('BLOCKED',path,','.join(issues))
    print('PUBLICATION_AUDIT', 'FAIL' if failures else 'PASS','snapshot_files',len(paths),'new_blobs',blobs)
    return bool(failures)
if __name__=='__main__':raise SystemExit(main())
