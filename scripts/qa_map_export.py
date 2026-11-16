"""Run synthetic map/order/sequencing QA."""
from synthetic_smoke import run
import tempfile
if __name__=='__main__':
    with tempfile.TemporaryDirectory(prefix='primer-map-qa-') as folder:run(folder)
