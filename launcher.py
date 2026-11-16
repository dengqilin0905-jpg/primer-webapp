"""Launch the installed application with its bundled public resources."""
import os,sys,subprocess
from pathlib import Path

def main():
    app=Path(__file__).resolve().parent/'app'
    env=dict(os.environ)
    env['PYTHONPATH']=str(app)+os.pathsep+env.get('PYTHONPATH','')
    return subprocess.call([sys.executable,'-m','streamlit','run',str(app/'app.py'),*sys.argv[1:]],env=env)
if __name__=='__main__':raise SystemExit(main())
