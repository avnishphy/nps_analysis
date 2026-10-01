(() => {
  'use strict';
  const PI = Math.PI, MP = 0.938, MPI0 = 0.1349768;
  const PP = [1.60077,-0.01523,37.08142,-4.11060,23.26192,0.00983,0.87073,-5.77115,-271.08678,0.13766,-0.00855,0.27885,-1.13212,-1.50415,-6.34766,0.55769,-0.01709];
  const PM = [1.75169,0.11144,47.35877,-4.69434,1.60552,0.00800,0.44194,-2.29188,-41.67194,0.69475,0.02527,-0.50178,-1.22825,-1.16878,5.75825,-1.00355,0.05055];
  const MODELS = {
    pi0_2021: {label:'2021 π⁰ proxy (active high-W core)', family:'2021', charge:'pi0', status:'ACTIVE CORE', validity:'Used by π⁰ branch at high W; blended with MAID below W=1.9 GeV.'},
    pip_2021: {label:'2021 π⁺ fit (active)', family:'2021', charge:'plus', status:'ACTIVE', validity:'Charged-pion fit; not a π⁰ model.'},
    pim_2021: {label:'2021 π⁻ fit (active)', family:'2021', charge:'minus', status:'ACTIVE', validity:'Charged-pion fit; not a π⁰ model.'},
    pip_blok: {label:'Blok π⁺ (dormant)', family:'blok', charge:'plus', status:'COMMENTED OUT', validity:'Historical charged-pion branch.'},
    pim_blok: {label:'Blok π⁻ (dormant)', family:'blok', charge:'minus', status:'COMMENTED OUT', validity:'Historical charged-pion branch.'},
    pip_04: {label:'param04 π⁺ (dormant)', family:'p04', charge:'plus', status:'COMMENTED OUT', validity:'Historical charged-pion branch; singular at |t|=0.02.'},
    pim_04: {label:'param04 π⁻ (dormant)', family:'p04', charge:'minus', status:'COMMENTED OUT', validity:'Historical charged-pion branch; singular at |t|=0.02.'},
    pip_3000: {label:'param_3000 π⁺ (dormant)', family:'p3000', charge:'plus', status:'COMMENTED OUT', validity:'Historical charged-pion branch.'},
    pim_3000: {label:'param_3000 π⁻ (dormant)', family:'p3000', charge:'minus', status:'COMMENTED OUT', validity:'Historical charged-pion branch.'}
  };
  const COLORS = ['#42d7bd','#ffc766','#55bfe0','#ff7b72','#c0a7ff','#80e080','#ff9fd3','#b7c8ff'];
  const $ = id => document.getElementById(id);
  const val = id => +$(id).value;
  const sq = x => x*x;
  function physical(Q2,W,u){
    const q0=(W*W-MP*MP-Q2)/(2*W), q=Math.sqrt(Math.max(0,q0*q0+Q2));
    const epi=(W*W+MPI0*MPI0-MP*MP)/(2*W), p=Math.sqrt(Math.max(0,epi*epi-MPI0*MPI0));
    const umin=Q2-MPI0*MPI0+2*q0*epi-2*q*p, umax=Q2-MPI0*MPI0+2*q0*epi+2*q*p;
    const c=q*p>0 ? (Q2-MPI0*MPI0+2*q0*epi-u)/(2*q*p) : 1;
    return {umin,umax,theta:Math.acos(Math.max(-1,Math.min(1,c))), valid:u>=umin-1e-9&&u<=umax+1e-9&&p>0};
  }
  function normParts(parts,k){ return {L:parts.L*k,T:parts.T*k,LT:parts.LT*k,TT:parts.TT*k}; }
  function fit2021(P,u,theta,Q2,S,pole){
    const fpi=pole/(1+P[0]*Q2+P[1]*Q2*Q2), A=Q2*fpi*fpi;
    return {
      L:(P[2]+P[14]/Q2)*u/sq(u+.02)*A*Math.exp(P[3]*u)/(S**P[10]+Math.sqrt(S)**P[16]),
      T:P[4]/Q2*Math.exp(P[5]*Q2*Q2)/(S**P[11]+Math.sqrt(S)**P[15])*Math.exp(P[13]*u),
      LT:P[6]/(1+P[9]*Q2)*Math.exp(P[7]*u)*Math.sin(theta)/S**P[12],
      TT:P[8]/(1+Q2)*Math.exp(-7*u)*sq(Math.sin(theta))
    };
  }
  function components(model,Q2,W,u,theta){
    const S=W*W, minus=model.charge==='minus', ch=minus?.25*(1+3*Math.exp(-10*u)):1;
    if(model.family==='2021'){
      if(model.charge==='pi0'){
        const a=fit2021(PP,u,theta,Q2,S,0), b=fit2021(PM,u,theta,Q2,S,0);
        const mean={}; for(const k of ['L','T','LT','TT']) mean[k]=(a[k]+b[k])/2;
        return normParts(mean,8.539/sq(S-MP*MP)/1e6);
      }
      return normParts(fit2021(model.charge==='plus'?PP:PM,u,theta,Q2,S,1),8.539/sq(S-MP*MP)/1e6);
    }
    let L,T,LT,TT,fpi;
    if(model.family==='blok'){
      L=27.8*Math.exp(-11.5*u); T=50*u*Math.exp(-5*u); LT=0; TT=-(4*L+.5*T)*sq(Math.sin(theta));
      fpi=1/(1+1.65*Q2+.5*Q2*Q2); L*=fpi*fpi*Q2/.1215; T=T/(.3+Q2)*ch; TT=TT/(.3+Q2)*ch;
      return normParts({L,T,LT,TT},15.333/sq(S-MP*MP)/1e6);
    }
    if(model.family==='p04'){
      fpi=1/(1+1.77*Q2+.05*Q2*Q2);
      L=350*Q2*fpi*fpi*Math.exp((16-7.5*Math.log(Q2))*(-u));
      T=(4.5/Q2+2/sq(Q2))*ch;
      LT=(Math.exp(.79-3.4/Math.sqrt(Q2)*u)+1.1-3.6/sq(Q2))*Math.sin(theta);
      TT=-5/sq(Q2)*(u/sq(u-.02))*sq(Math.sin(theta))*ch;
    } else {
      fpi=1/(1+1.77*Q2+.12*Q2*Q2);
      L=19.8*u/sq(u+.02)*Q2*fpi*fpi*Math.exp(-3.66*u);
      T=5.96/Q2*Math.exp(-.013*Q2*Q2)*ch;
      LT=16.533/(1+Q2)*Math.exp(-5.1437*u)*Math.sin(theta);
      TT=-178.06/(1+Q2)*Math.exp(-7.1381*u)*sq(Math.sin(theta))*ch;
    }
    return normParts({L,T,LT,TT},8.539/sq(S-MP*MP)/1e6);
  }
  function evaluate(model,Q2,W,u,theta,phi,eps){
    const c=components(model,Q2,W,u,theta), a=Math.sqrt(2*eps*(1+eps));
    const terms={U:c.T+eps*c.L,LT:a*Math.cos(phi)*c.LT,TT:eps*Math.cos(2*phi)*c.TT};
    return {...c,...terms,total:(terms.U+terms.LT+terms.TT)/(2*PI)};
  }
  const toNb=x=>x*1e9;
  function fmt(x){ if(!Number.isFinite(x)) return String(x); const a=Math.abs(x); return a!==0&&(a<.001||a>=1e5)?x.toExponential(4):x.toFixed(a<1?5:a<100?3:1); }
  function state(overrides={}){
    const model=MODELS[$('model').value], Q2=overrides.Q2??val('q2'), W=overrides.W??val('w'), eps=overrides.eps??val('eps');
    const kin0=physical(Q2,W,0), tp=overrides.tp??val('tp'), u=kin0.umin+tp, kin=physical(Q2,W,u);
    const theta=$('physical-angle').checked?kin.theta:val('theta')*PI/180, phi=(overrides.phiDeg??val('phi'))*PI/180;
    return {model,Q2,W,eps,tp,u,theta,phi,kin};
  }
  function check(s){
    const problems=[];
    if(!(s.Q2>0)) problems.push('Q² must be positive.'); if(!(s.W>MP+MPI0)) problems.push('W is below the π⁰p threshold.');
    if(!(s.eps>=0&&s.eps<=1)) problems.push('ε must lie from 0 to 1.'); if(!s.kin.valid&&$('physical-angle').checked) problems.push('−t′ is outside the physical two-body range.');
    return problems;
  }
  function formula(m){
    if(m.family==='2021') return `2021 response fit (physics_pion.f: sig_param_2021)\nF = [T + εL + √(2ε(1+ε)) LT cosφ + ε TT cos(2φ)]/(2π)\nnormalization = 8.539 / (W² − 0.938²)² / 10⁶\n${m.charge==='pi0'?'π⁰ = arithmetic mean of π⁺ and π⁻ calls; pole factor fπ forced to 0.':'charged fit; pole factor fπ enabled.'}`;
    if(m.family==='blok') return `Dormant sig_blok branch\nL₀=27.8e^(−11.5|t|), T₀=50|t|e^(−5|t|), LT=0\nTT₀=−(4L₀+0.5T₀)sin²θ; fπ=1/(1+1.65Q²+0.5Q⁴)\nnormalization = 15.333 / (W² − M²)² / (2π·10⁶)`;
    if(m.family==='p04') return `Dormant sig_param04 branch (SIMC passes positive t=|t|)\nL=350Q²fπ² exp[(16−7.5lnQ²)(−t)]\nT=4.5/Q²+2/Q⁴; LT=[e^(0.79−3.4t/√Q²)+1.1−3.6/Q⁴]sinθ\nTT=−5/Q⁴ · t/(t−0.02)² · sin²θ`;
    return `Dormant sig_param_3000 branch\nL=19.8t/(t+0.02)² Q²fπ²e^(−3.66t); T=5.96/Q² e^(−0.013Q⁴)\nLT=16.533/(1+Q²)e^(−5.1437t)sinθ\nTT=−178.06/(1+Q²)e^(−7.1381t)sin²θ`;
  }
  function linesForPhi(s){
    const xs=Array.from({length:181},(_,i)=>i*2), series=[{name:'total',color:COLORS[0],y:[]},{name:'T + εL',color:COLORS[1],y:[]},{name:'LT modulation',color:COLORS[2],y:[]},{name:'TT modulation',color:COLORS[3],y:[]}];
    for(const x of xs){const r=evaluate(s.model,s.Q2,s.W,s.u,s.theta,x*PI/180,s.eps); series[0].y.push(toNb(r.total)); series[1].y.push(toNb(r.U/(2*PI))); series[2].y.push(toNb(r.LT/(2*PI))); series[3].y.push(toNb(r.TT/(2*PI)));}
    return {xs,series,xlabel:'φ (degrees)',ylabel:'d²σ/dt dφ  [nb/GeV²/rad]'};
  }
  function linesForT(s){
    const max=Math.max(.05,Math.min(1.5,s.kin.umax-s.kin.umin)), xs=Array.from({length:121},(_,i)=>max*i/120), keys=[['total','d²σ at selected φ'],['U','σU = σT + εσL'],['LT','σLT'],['TT','σTT']];
    const series=keys.map((k,i)=>({name:k[1],key:k[0],color:COLORS[i],y:[]}));
    for(const x of xs){const k=physical(s.Q2,s.W,s.kin.umin+x), r=evaluate(s.model,s.Q2,s.W,s.kin.umin+x,k.theta,s.phi,s.eps); for(const a of series)a.y.push(toNb(a.key==='total'?r.total:r[a.key]));}
    return {xs,series,xlabel:'−t′ = |t| − |t|min  [GeV²]',ylabel:'response [nb/GeV²]  (total is /rad)'};
  }
  function linesForQ(s){
    const xs=Array.from({length:121},(_,i)=>.3+5.7*i/120), series=[{name:'d²σ at selected φ',color:COLORS[0],y:[]},{name:'σU',color:COLORS[1],y:[]},{name:'σLT',color:COLORS[2],y:[]},{name:'σTT',color:COLORS[3],y:[]}];
    for(const x of xs){const k0=physical(x,s.W,0), u=k0.umin+s.tp, k=physical(x,s.W,u); if(!k.valid){series.forEach(a=>a.y.push(NaN));continue;} const r=evaluate(s.model,x,s.W,u,k.theta,s.phi,s.eps); [r.total,r.U,r.LT,r.TT].forEach((v,i)=>series[i].y.push(toNb(v)));}
    return {xs,series,xlabel:'Q² [GeV²]',ylabel:'response [nb/GeV²]  (total is /rad)'};
  }
  function nice(v){if(!Number.isFinite(v)||v===0)return 1; const p=10**Math.floor(Math.log10(Math.abs(v))), n=Math.abs(v)/p; return (n<=1?1:n<=2?2:n<=5?5:10)*p;}
  function draw(id,data){
    const host=$(id), W=820,H=330,m={l:84,r:20,t:22,b:56}; let ys=[]; data.series.forEach(s=>ys.push(...s.y.filter(Number.isFinite))); let ymin=Math.min(...ys), ymax=Math.max(...ys); if(!Number.isFinite(ymin)||ymin===ymax){ymin=0;ymax=1;} const pad=(ymax-ymin)*.08; ymin-=pad;ymax+=pad; if(ymin>0)ymin=0; const x0=data.xs[0],x1=data.xs.at(-1), X=x=>m.l+(x-x0)/(x1-x0)*(W-m.l-m.r),Y=y=>m.t+(ymax-y)/(ymax-ymin)*(H-m.t-m.b);
    let svg=`<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${data.xlabel} plot">`;
    for(let i=0;i<=5;i++){const y=ymin+(ymax-ymin)*i/5,py=Y(y);svg+=`<line class="gridline" x1="${m.l}" y1="${py}" x2="${W-m.r}" y2="${py}"/><text x="${m.l-9}" y="${py+4}" text-anchor="end">${fmt(y)}</text>`;}
    for(let i=0;i<=6;i++){const x=x0+(x1-x0)*i/6,px=X(x);svg+=`<line class="gridline" x1="${px}" y1="${m.t}" x2="${px}" y2="${H-m.b}"/><text x="${px}" y="${H-m.b+19}" text-anchor="middle">${fmt(x)}</text>`;}
    svg+=`<line class="axis" x1="${m.l}" y1="${H-m.b}" x2="${W-m.r}" y2="${H-m.b}"/><line class="axis" x1="${m.l}" y1="${m.t}" x2="${m.l}" y2="${H-m.b}"/><text class="axis-title" x="${(m.l+W-m.r)/2}" y="${H-9}" text-anchor="middle">${data.xlabel}</text><text class="axis-title" transform="translate(16 ${(m.t+H-m.b)/2}) rotate(-90)" text-anchor="middle">${data.ylabel}</text>`;
    for(const s of data.series){let d='',pen=false;s.y.forEach((y,i)=>{if(!Number.isFinite(y)){pen=false;return;}d+=(pen?'L':'M')+X(data.xs[i]).toFixed(2)+' '+Y(y).toFixed(2)+' ';pen=true;});svg+=`<path class="trace" stroke="${s.color}" d="${d}"/>`;}
    svg+=`<rect class="hit" x="${m.l}" y="${m.t}" width="${W-m.l-m.r}" height="${H-m.t-m.b}" data-hit="${id}"/></svg><div class="legend">${data.series.map(s=>`<span><i style="background:${s.color}"></i>${s.name}</span>`).join('')}</div>`; host.innerHTML=svg; host._plot={data,W,H,m};
    const hit=host.querySelector('.hit'),tip=host.closest('.chart-card').querySelector('.plot-tooltip'); hit.addEventListener('mousemove',e=>{const r=hit.getBoundingClientRect(), px=(e.clientX-r.left)/r.width, x=x0+px*(x1-x0), idx=Math.max(0,Math.min(data.xs.length-1,Math.round(px*(data.xs.length-1))));tip.innerHTML=`<b>${data.xlabel.split('[')[0]} ${fmt(data.xs[idx])}</b><br>`+data.series.map(s=>`${s.name}: ${fmt(s.y[idx])}`).join('<br>');tip.style.display='block';tip.style.left=Math.min(e.offsetX+18,host.clientWidth-250)+'px';tip.style.top=(e.offsetY+55)+'px';}); hit.addEventListener('mouseleave',()=>tip.style.display='none');
  }
  function csv(data){let out=[['x',...data.series.map(s=>s.name)].join(',')];data.xs.forEach((x,i)=>out.push([x,...data.series.map(s=>s.y[i])].join(',')));return out.join('\n');}
  function download(name,text){const a=document.createElement('a');a.href=URL.createObjectURL(new Blob([text],{type:'text/csv'}));a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(a.href),500);}
  function update(){
    const s=state(), problems=check(s), r=evaluate(s.model,s.Q2,s.W,s.u,s.theta,s.phi,s.eps), p=[linesForPhi(s),linesForT(s),linesForQ(s)];
    $('error').textContent=problems.join(' '); $('error').classList.toggle('show',problems.length>0); $('model-status').textContent=s.model.status; $('model-status').className='status-chip '+(s.model.status.startsWith('ACTIVE')?'good':'warn'); $('model-note').textContent=s.model.validity;
    $('kin-status').textContent=s.kin.valid?'physical 2-body point':'outside physical range'; $('kin-status').className='status-chip '+(s.kin.valid?'good':'warn');
    $('theta').disabled=$('physical-angle').checked; $('theta').value=(s.theta*180/PI).toFixed(3); $('theta-readout').textContent=(s.theta*180/PI).toFixed(3); $('umin').textContent=fmt(s.kin.umin); $('umax').textContent=fmt(s.kin.umax); $('uvalue').textContent=fmt(s.u);
    $('m-total').textContent=fmt(toNb(r.total)); $('m-u').textContent=fmt(toNb(r.U)); $('m-lt').textContent=fmt(toNb(r.LT)); $('m-tt').textContent=fmt(toNb(r.TT)); $('formula').textContent=formula(s.model);
    draw('phi-chart',p[0]);draw('t-chart',p[1]);draw('q-chart',p[2]); $('phi-chart')._csv=p[0];$('t-chart')._csv=p[1];$('q-chart')._csv=p[2];
  }
  document.querySelectorAll('[data-control]').forEach(x=>x.addEventListener('input',update)); document.querySelectorAll('[data-csv]').forEach(b=>b.addEventListener('click',()=>download(b.dataset.csv+'.csv',csv($(b.dataset.target)._csv))));
  update();
})();
