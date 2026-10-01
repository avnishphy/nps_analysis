"""Archive public source candidates only; never open analysis ROOT files."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
import urllib.request

BASE = Path('/tmp/nps_livetime_4398_sources_20260909')
BASE.mkdir(exist_ok=True)
SOURCES = {
    'doc1001.html': 'https://hallcweb.jlab.org/doc-public/ShowDocument?docid=1001',
    'doc1022.html': 'https://hallcweb.jlab.org/doc-public/ShowDocument?docid=1022',
    'doc1028.html': 'https://hallcweb.jlab.org/doc-public/ShowDocument?docid=1028',
    'doc1063.html': 'https://hallcweb.jlab.org/doc-public/ShowDocument?docid=1063',
    'doc1110.html': 'https://hallcweb.jlab.org/doc-public/ShowDocument?docid=1110',
    'nps_trigger_status_20230202.pdf': 'https://wiki.jlab.org/cuawiki/images/9/90/NPS_Trigger_Status_Feb_2_2023.pdf',
    'nps_pi0_status_20250506.pdf': 'https://indico.jlab.org/event/946/contributions/16518/attachments/12601/20084/Edited_SemiInclusive_Pi0_Status_May2025_styled_clean.pdf',
    'shms_instrument_draft_20240605.pdf': 'https://mailman.jlab.org/pipermail/hallc/attachments/20240605/2611cf32/attachment-0001.pdf',
    'f2_livetime_20220217.pdf': 'https://indico.jlab.org/event/517/contributions/9310/attachments/7536/10488/F2_HallCwinterCollaboration_2022.pdf',
    'hallc_daq.html': 'https://hallcweb.jlab.org/wiki/index.php?title=Hall_C_DAQ',
    'hallc_trigger_layout.html': 'https://hallcweb.jlab.org/wiki/index.php?title=Hall_C_Trigger_Layout',
    'hallc_scalers.html': 'https://hallcweb.jlab.org/wiki/index.php?title=Scalers',
    'hallc_livetime_studies.html': 'https://hallcweb.jlab.org/wiki/index.php?title=Electronic%2FComputer_Live_Time_Studies',
    'nps_analysis_index.html': 'https://wiki.jlab.org/cuawiki/index.php/NPS_RG1a_Analysis',
    'nps_replay_readme.md': 'https://raw.githubusercontent.com/JeffersonLab/nps_replay/develop/README.md',
}

def fetch(item):
    name, url = item
    record = dict(file=name, url=url)
    try:
        with urllib.request.urlopen(url, timeout=35) as response:
            data = response.read()
            record.update(final_url=response.url, content_type=response.headers.get('Content-Type'))
        if name.endswith('.pdf') and not data.startswith(b'%PDF-'):
            raise ValueError('Response is not a PDF')
        path = BASE / name
        if path.exists() and path.read_bytes() != data:
            raise ValueError('Refusing to overwrite existing different content')
        path.write_bytes(data)
        record.update(status='downloaded', bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
    except Exception as exc:
        record.update(status='failed', error=str(exc))
    return record

if __name__ == '__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        records = list(pool.map(fetch, SOURCES.items()))
    (BASE / 'download_manifest.json').write_text(json.dumps(records, indent=2) + '\n')
    for record in records:
        print(record['file'], record['status'], record.get('bytes', record.get('error')))
