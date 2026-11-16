"""Optional private resources. Importing software never creates or discovers lab data."""
import json,os
from pathlib import Path
REPO_ROOT=Path(__file__).resolve().parents[1]

def private_settings():
    filename=os.environ.get('PRIMER_PRIVATE_CONFIG','').strip()
    if not filename:return {}
    path=Path(filename).expanduser().resolve()
    with path.open(encoding='utf-8') as handle:data=json.load(handle)
    if not isinstance(data,dict):raise ValueError('Private resource configuration must be a JSON object')
    return {key:str((path.parent/Path(value).expanduser()).resolve()) for key,value in data.items() if value not in (None,'')}

def resource_path(key,default):
    value=os.environ.get('PRIMER_'+key.upper()) or private_settings().get(key)
    return Path(value).expanduser().resolve() if value else Path(default)

def state_root():
    return resource_path('state_dir',Path.home()/'.local/share/primer-webapp')
