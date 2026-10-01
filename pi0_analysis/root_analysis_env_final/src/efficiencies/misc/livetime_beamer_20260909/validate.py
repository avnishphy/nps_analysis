"""Validate slide layout, embedded numerical tables and frozen evidence integrity."""
from pathlib import Path
import subprocess,xml.etree.ElementTree as ET,json,hashlib,re,math
from PIL import Image,ImageDraw
P=Path(__file__).resolve().parent;B=P/'build';pdf=P/'KinC_x60_4b_LH2_livetime_Beamer.pdf'
subprocess.run(['pdftotext','-bbox',str(pdf),str(B/'bbox.html')],check=True)
raw=(B/'bbox.html').read_text()
root=ET.fromstring(re.sub(r'[\x00-\x08\x0b\x0c\x0e-\x1f]','',raw));ns={'x':'http://www.w3.org/1999/xhtml'}
pages=root.findall('.//x:page',ns);issues=[];titles=[]
for i,p in enumerate(pages,1):
    w,h=float(p.attrib['width']),float(p.attrib['height'])
    words=p.findall('x:word',ns)
    if not words:issues.append({'page':i,'error':'no text'})
    for t in words:
        a=t.attrib
        if float(a['xMin']) < -0.5 or float(a['yMin']) < -0.5 or float(a['xMax']) > w+.5 or float(a['yMax']) > h+.5:
            issues.append({'page':i,'error':'out of page','text':t.text,'box':a})
    title=' '.join(t.text or '' for t in words if float(t.attrib['yMin'])<40)
    titles.append({'page':i,'title':title})
log=(B/'livetime.log').read_text();layout_errors=[s for s in log.splitlines() if 'Overfull' in s or s.startswith('!')]
R=json.loads((P/'baseline/audited_results.json').read_text());prod=[r for r in R if r['status']=='CALCULATED' and r['run_type']=='production']
for r in prod:
    n=r['nominal'];p=r['prescale_factor']
    for a,b in [(n['EDTM_tight'],p*n['E_tight']/n['D']),(n['CLT_physics'],p*(n['N']-n['E_tight'])/(n['S']-n['D'])),(n['CLT_all'],(1-n['D']/n['S'])*n['CLT_physics']+n['D']/n['S']*n['EDTM_tight'])]:
        assert math.isclose(a,b,rel_tol=0,abs_tol=1e-12),(r['run'],a,b)
tables='\n'.join(t.read_text() for t in (P/'tables').glob('*.tex'))
for r in prod:
    assert str(r['run']) in tables
    for v in [r['old_ratio'],r['nominal']['EDTM_raw'],r['nominal']['EDTM_tight'],r['nominal']['CLT_physics']]:assert f'{v:.6f}' in tables
manifest=json.loads((P/'sources/local_manifest.json').read_text())
for r in manifest:assert hashlib.sha256((P/r['file']).read_bytes()).hexdigest()==r['sha256'],r['file']
for r in json.loads((P/'sources/web_manifest.json').read_text()):
    if r['status']=='downloaded':assert hashlib.sha256((P/'sources'/r['file']).read_bytes()).hexdigest()==r['sha256'],r['file']
tex=(P/'livetime.tex').read_text()
required=['corrected-time tails','No PRE-derived correction','same-cache','good-helicity','time-weighted','diagnostic']
text=(B/'slides.txt').read_text().lower()
for term in required:assert term.lower() in text,term
first_backup=next(t['page'] for t in titles if t['title'].startswith('Backup:'))
result={'pages':len(pages),'main_slides':first_backup-1,'backup_slides':len(pages)-first_backup+1,'out_of_page_text':issues,'overfull_or_tex_errors':layout_errors,'production_formula_checks':len(prod),'production_table_rows':len(prod),'copied_evidence_hashes_checked':len(manifest),'all_source_hashes_checked':True,'titles':titles}
(P/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
assert not issues and not layout_errors,result
preview=B/'previews';preview.mkdir(exist_ok=True)
for start in range(1,len(pages)+1,12):
    im=Image.new('RGB',(1280,4*385),'#d9e1e6');draw=ImageDraw.Draw(im)
    for k,i in enumerate(range(start,min(start+12,len(pages)+1))):
        f=preview/f'slide-{i:02d}.png'
        if not f.exists():continue
        pic=Image.open(f).convert('RGB');pic.thumbnail((426,355))
        x=(k%3)*426;y=(k//3)*385
        im.paste(pic,(x,y+22));draw.text((x+8,y+4),f'Page {i}',fill='#14334a')
    im.save(preview/f'contact_{start:02d}.png')
print(json.dumps({k:v for k,v in result.items() if k!='titles'},indent=2))
