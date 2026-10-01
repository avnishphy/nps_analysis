"""Fetch public PDF attachments linked from the archived DocDB metadata."""
from nps_livetime_collect import BASE, fetch
import concurrent.futures
import html
import json
import re
from urllib.parse import urlsplit, unquote

items = []
for path in sorted(BASE.glob('doc*.html')):
    for raw in re.findall(r'href="([^"]+)"', path.read_text()):
        url = html.unescape(raw)
        if url.startswith('https://hallcweb.jlab.org/DocDB/') and url.endswith('.pdf'):
            items.append((path.stem + '_' + unquote(urlsplit(url).path.rsplit('/', 1)[-1]).replace(' ', '_'), url))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    records = list(pool.map(fetch, items))
(BASE / 'attachment_manifest.json').write_text(json.dumps(records, indent=2) + '\n')
for record in records:
    print(record['file'], record['status'], record.get('bytes', record.get('error')))
