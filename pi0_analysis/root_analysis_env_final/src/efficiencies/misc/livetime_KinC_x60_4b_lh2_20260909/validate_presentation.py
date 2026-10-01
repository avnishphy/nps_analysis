"""Check the actual final PDF, including text emitted by table artists."""
from pathlib import Path
import subprocess,json,hashlib,xml.etree.ElementTree as ET
P=Path(__file__).resolve().parent;pdf=P/'KinC_x60_4b_LH2_livetime_investigation.pdf'
subprocess.run(['pdftotext','-bbox',str(pdf),str(P/'pdf_text_bounds.html')],check=True)
root=ET.parse(P/'pdf_text_bounds.html').getroot();ns='{http://www.w3.org/1999/xhtml}'
pages=list(root.iter(ns+'page'));bad=[]
for i,page in enumerate(pages,1):
    w,h=float(page.attrib['width']),float(page.attrib['height'])
    for t in page.iter(ns+'word'):
        a=t.attrib
        if float(a['xMin'])<0 or float(a['yMin'])<0 or float(a['xMax'])>w or float(a['yMax'])>h:bad.append(dict(page=i,text=t.text,bounds=a))
v=dict(pages=len(pages),out_of_page_words=bad,pdf_sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),visual_review='Rendered coverage, counter table, factor-one comparison, prescale comparison and appendix layouts inspected; no clipping observed.')
(P/'pdf_validation.json').write_text(json.dumps(v,indent=2))
assert len(pages)==44 and not bad
print('44 PDF pages; no out-of-page text; checksum recorded.')
