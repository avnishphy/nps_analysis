"""Archive public CODA primary sources; no login, staging or external messages."""
from pathlib import Path
import urllib.request,json,hashlib,datetime,concurrent.futures,sys
P=Path(__file__).resolve().parent;S=P/'sources';S.mkdir(exist_ok=True)
items={
 'ti_index.html':'https://coda.jlab.org/drupal/system/files/coda/LibraryManual/tiLib/index.html',
 'ti_files.html':'https://coda.jlab.org/drupal/system/files/coda/LibraryManual/tiLib/files.html',
 'tiLib_8c.html':'https://coda.jlab.org/drupal/system/files/coda/LibraryManual/tiLib/tiLib_8c.html',
 'TI.pdf':'https://coda.jlab.org/drupal/system/files/pdfs/HardwareManual/TI/TI.pdf',
 'TS.pdf':'https://coda.jlab.org/drupal/system/files/pdfs/HardwareManual/TS/TS.pdf',
 'single_crate.html':'https://coda.jlab.org/drupal/content/example-single-crate-configuration',
 'walkthrough.html':'https://coda.jlab.org/drupal/system/files/coda/3.10/walkthrough/index.html',
 'roc.html':'https://coda.jlab.org/drupal/content/readout-controller-roc',
 'documents.html':'https://coda.jlab.org/drupal/content/documents',
 'vme_drivers.html':'https://coda.jlab.org/drupal/content/vme-module-drivers',
}
def fetch(pair):
 name,url=pair;row=dict(file=name,url=url,utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
 try:
  with urllib.request.urlopen(url,timeout=40) as r:b=r.read();row.update(status=r.status,final_url=r.url,content_type=r.headers.get('Content-Type'))
  (S/name).write_bytes(b);row.update(bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
 except Exception as e:row['error']=str(e)
 print(name,row.get('status'),row.get('bytes'),row.get('error',''),flush=True);return row
manifest=P/'source_manifest.json'
previous=json.loads(manifest.read_text()) if manifest.exists() else []
if len(sys.argv)>1:
 items={url.rsplit('/',1)[-1]:url for url in sys.argv[1:]}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:rows=list(pool.map(fetch,items.items()))
replaced={r['file'] for r in rows}
manifest.write_text(json.dumps([r for r in previous if r['file'] not in replaced]+rows,indent=2)+'\n')
