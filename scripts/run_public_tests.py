"""Portable synthetic regression runner; never connects to a laboratory database."""
import os,sys,tempfile,unittest
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(root/'app'),str(root/'test')]
with tempfile.TemporaryDirectory(prefix='primer-public-') as folder:
    for key in list(os.environ):
        if key.startswith('PRIMER_') or key in ('DB_PATH','VERIFICATION_STATE_DIR'):os.environ.pop(key)
    os.environ['HOME']=folder
    os.environ['DB_PATH']=str(Path(folder)/'app.db')
    os.environ['VERIFICATION_STATE_DIR']=str(Path(folder)/'verification')
    os.environ['PYTHONDONTWRITEBYTECODE']='1'
    suite=unittest.defaultTestLoader.discover(str(root/'test'))
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    if result.skipped:raise SystemExit('Skipped public tests are not acceptance evidence')
    raise SystemExit(not result.wasSuccessful())
