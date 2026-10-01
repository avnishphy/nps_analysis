"""Read public sources; record bytes, URLs and retrieval failures. No login/staging."""
from pathlib import Path
import urllib.request, hashlib, json, datetime, subprocess
P=Path(__file__).resolve().parent/'sources'
URLS={
 'zhang_deadtime_20240219.pdf':'https://redmine.jlab.org/attachments/download/2343/Dead%20time%20analysis_2024-02-19.pdf',
 'zhang_deadtime_20240303.pdf':'https://redmine.jlab.org/attachments/download/2368/Updates_2024-03-03.pdf',
 'zhang_deadtime_20240718.pdf':'https://indico.jlab.org/event/866/contributions/14918/subcontributions/283/attachments/11529/17843/Deadtime%20and%20Efficiency.pdf',
 'raydo_trigger_20240717.pdf':'https://indico.jlab.org/event/866/contributions/14914/subcontributions/234/attachments/11505/17805/NPS_July17_2024_Raydo.pdf',
 'pionlt_EDTM_study.pdf':'https://redmine.jlab.org/attachments/download/1499/EDTM_Study_Report.pdf',
 'issue866.html':'https://redmine.jlab.org/issues/866',
 'issue837.html':'https://redmine.jlab.org/issues/837',
 'issue836.html':'https://redmine.jlab.org/issues/836',
 'trigger_history.html':'https://hallcweb.jlab.org/wiki/index.php/Trigger_History',
 'edtm_pulser.html':'https://hallcweb.jlab.org/wiki/index.php?title=Hall_C_EDTM_Pulser',
}
rows=json.loads((P/'web_manifest.json').read_text()) if (P/'web_manifest.json').exists() else []
for name,url in URLS.items():
    if any(r['file']==name and r['status']=='downloaded' for r in rows):continue
    row={'file':name,'url':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    try:
        with urllib.request.urlopen(url,timeout=25) as response:
            data=response.read();row['final_url']=response.url
        if name.endswith('.pdf') and not data.startswith(b'%PDF'):raise ValueError('Not a PDF')
        (P/name).write_bytes(data)
        row.update(bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),status='downloaded')
        if name.endswith('.pdf'):subprocess.run(['pdftotext','-layout',str(P/name),str(P/Path(name).with_suffix('.txt'))],check=True)
    except Exception as exc:row.update(status='failed',error=str(exc))
    rows.append(row);print(name,row['status'],row.get('bytes',row.get('error')),flush=True)
    (P/'web_manifest.json').write_text(json.dumps(rows,indent=2)+'\n')
