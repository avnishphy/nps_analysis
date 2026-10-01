from nps_livetime_collect import BASE, fetch
import json
items = [
 ('singh_nps_luminosity_202505.pdf', 'https://indico.jlab.org/event/946/contributions/16510/attachments/12602/20074/Luminosity_nps_collaboration_meeting_2025.pdf'),
 ('yero_EDTM_report_20170906.pdf', 'https://raw.githubusercontent.com/JeffersonLab/Hall-C-Trigger-Setup/master/edtm_studies/EDTM_report.pdf'),
 ('yero_EDTM_LT_Studies_20171121.pdf', 'https://raw.githubusercontent.com/JeffersonLab/Hall-C-Trigger-Setup/master/edtm_studies/EDTM_LT_Studies.pdf'),
 ('nps_runplan_archive_20230911.pdf', 'https://hallcweb.jlab.org/wiki/images/archive/1/10/20230911123553%21NPS_DVCS_RunPlan.pdf'),
]
records = [fetch(item) for item in items]
(BASE / 'supplement_manifest.json').write_text(json.dumps(records, indent=2) + '\n')
for r in records:
    print(r['file'], r['status'], r.get('bytes', r.get('error')))
