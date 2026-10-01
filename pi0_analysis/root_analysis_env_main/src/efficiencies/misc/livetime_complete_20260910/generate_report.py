"""Generate an editable article body and tables. Not run by build.sh."""
from pathlib import Path
import csv,re
P=Path(__file__).resolve().parent
def esc(s):
 return ''.join({'\\':r'\textbackslash{}','&':r'\&','%':r'\%','$':r'\$','#':r'\#','_':r'\_','{':r'\{','}':r'\}','~':r'\textasciitilde{}','^':r'\textasciicircum{}'}.get(c,c) for c in s)
def inline(s):
 pat=r'(\$[^$]+\$|`[^`]+`|\[[^\]]+\]\([^)]+\)|\*\*[^*]+\*\*)';out=[]
 for i,t in enumerate(re.split(pat,s)):
  if not i%2:out.append(esc(t));continue
  if t.startswith('$'):out.append(t)
  elif t.startswith('`'):out.append(r'\path|'+t[1:-1]+'|')
  elif t.startswith('**'):out.append(r'\textbf{'+esc(t[2:-2])+'}')
  else:
   m=re.fullmatch(r'\[([^]]+)\]\(([^)]+)\)',t);url=m[2].replace('%',r'\%')
   out.append(r'\href{'+url+'}{'+esc(m[1])+'}')
 return ''.join(out)
def render(text):
 lines=text.splitlines();out=[];i=0
 while i<len(lines):
  line=lines[i].strip()
  if not line:i+=1;continue
  if line=='$$':
   j=i+1
   while lines[j].strip()!='$$':j+=1
   out.append('\\[\n'+'\n'.join(lines[i+1:j])+'\n\\]');i=j+1;continue
  if line.startswith('```'):
   j=i+1
   while j<len(lines) and not lines[j].startswith('```'):j+=1
   out.append('\\begin{lstlisting}\n'+'\n'.join(lines[i+1:j])+'\n\\end{lstlisting}');i=j+1;continue
  if line.startswith('#'):
   n=len(line)-len(line.lstrip('#'));title=re.sub(r'^\d+\.\s+', '', line[n:].strip())
   if n>1:out.append('\\'+('section' if n==2 else 'subsection')+'{'+inline(title)+'}')
   i+=1;continue
  if line.startswith('!['):
   m=re.match(r'!\[([^]]*)\]\(([^)]+)\)',line);f=m[2]
   if f.endswith('.png') and (P/f[:-4]).with_suffix('.pdf').exists():f=f[:-4]+'.pdf'
   out.append(r'\begin{figure}[htbp]\centering\includegraphics[width=\linewidth,height=.72\textheight,keepaspectratio]{'+f+r'}\caption{'+inline(m[1])+r'}\end{figure}\clearpage');i+=1;continue
  if line.startswith('|'):
   rows=[]
   while i<len(lines) and lines[i].strip().startswith('|'):
    cells=[c.strip() for c in lines[i].strip().strip('|').split('|')]
    if not all(re.fullmatch(r'[-: ]+',c) for c in cells):rows.append(cells)
    i+=1
   n=len(rows[0]);width=(1-.014*(n-1))/n
   if n==3:widths=[.12,.43,.422] if rows[0][0]=='ID' else [.24,.34,.392]
   elif n==4:widths=[.18,.26,.24,.278]
   else:widths=[int(width*1000)/1000]*n
   spec='@{}'+'@{\\hspace{.014\\linewidth}}'.join('p{'+f'{w:.3f}'+r'\linewidth}' for w in widths)+'@{}'
   out.append('{\\small\\setlength{\\tabcolsep}{0pt}\\renewcommand{\\arraystretch}{1.15}\n\\begin{longtable}{'+spec+'}\\toprule\n'+' & '.join(r'\textbf{'+inline(c)+'}' for c in rows[0])+r'\\\midrule\endhead')
   for row in rows[1:]:out.append(' & '.join(inline(c) for c in row)+r'\\')
   out.append('\\bottomrule\\end{longtable}}');continue
  if re.match(r'^(- |\d+\. )',line):
   ordered=bool(re.match(r'^\d+\.',line));env='enumerate' if ordered else 'itemize';out.append('\\begin{'+env+'}')
   while i<len(lines) and re.match(r'^(- |\d+\. )',lines[i].strip()):
    item=re.sub(r'^(- |\d+\. )','',lines[i].strip());i+=1
    while i<len(lines) and lines[i].strip() and not re.match(r'^(- |\d+\. |#|\|)',lines[i].strip()):item+=' '+lines[i].strip();i+=1
    out.append('\\item '+inline(item))
   out.append('\\end{'+env+'}');continue
  para=line;i+=1
  while i<len(lines) and lines[i].strip() and not re.match(r'^(#|\||!\[|\$\$|```|- |\d+\. )',lines[i].strip()):para+=' '+lines[i].strip();i+=1
  out.append(inline(para)+'\n\n')
 return '\n'.join(out)
(P/'sections/report_body.tex').write_text(render((P/'REPORT.md').read_text()))
(P/'sections/references.tex').write_text(render((P/'REFERENCES.md').read_text()))
rows=list(csv.DictReader((P/'data/production_summary.csv').open()))
count=['\\section{All 43 production count inputs}',r'\small\begin{longtable}{rrrrrrrr}\toprule Run & TI & p & Seg. & N & Eraw & D & S\\\midrule\endhead']
rat=['\\section{All 43 production ratio comparisons}',r'\small\begin{longtable}{rrrrrrr}\toprule Run & Same-file old & Matched raw & Tight & pN/S & Proposed phys. & Flag\\\midrule\endhead']
for r in rows:
 count.append(' & '.join(r[k] if k not in ['D','S'] else str(int(float(r[k]))) for k in ['run','trigger','p','segments','N','E_raw','D','S'])+r'\\')
 rat.append(r['run']+' & '+' & '.join(f'{float(r[k]):.7f}' for k in ['same_cache_NewGen','matched_raw_EDTM','matched_tight_EDTM','all_trigger_ratio','raw_subtracted_physics_ratio'])+' & '+('!' if int(r['counter_suspect_intervals']) else '--')+r'\\')
for a in [count,rat]:a.append(r'\bottomrule\end{longtable}\normalsize')
(P/'sections/tables.tex').write_text('\n'.join(count+['\\clearpage']+rat))
for j in range(0,len(rows),9):
 s=[r'\begin{tabular}{rrrrrrr}\toprule Run & TI/p & Seg. & Old & Raw & Tight & Flag\\\midrule']
 for r in rows[j:j+9]:s.append(r['run']+' & '+r['trigger']+'/'+r['p']+' & '+r['segments']+' & '+' & '.join(f'{float(r[k]):.6f}' for k in ['same_cache_NewGen','matched_raw_EDTM','matched_tight_EDTM'])+' & '+('!' if int(r['counter_suspect_intervals']) else '--')+r'\\')
 s.append(r'\bottomrule\end{tabular}');(P/'sections'/f'slide_table_{j//9+1}.tex').write_text('\n'.join(s))
print('Generated editable report body, reference appendix and 43-run tables.')
