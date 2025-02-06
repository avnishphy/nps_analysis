 // Filename: PeterB.C

// Makes skim files with needed variables for good el. in HMS
// and pi0 in NPS. Also makes a few histograms
// Needs runnumber
#include <stdio.h>  

void PeterB(Int_t runNumber){

// Set up the Skim files
  TString fileNameSkim = 
     "/w/hallc-scshelf2102/nps/bosted/NewSkimfiles/Skim" ;
  TString fileNameChist = 
    "/w/hallc-scshelf2102/nps/bosted/NewSkimfiles/Hhist" ;
  TString fileNameSkimnps = 
    "/w/hallc-scshelf2102/nps/bosted/NewSkimfiles/Skimnps" ;

  fileNameSkim += runNumber; 
  fileNameChist += runNumber; 
  fileNameSkimnps += runNumber; 

  fileNameSkim += ".txt"; 
  fileNameChist += ".txt";
  fileNameSkimnps += ".txt"; 

  FILE *f2 = fopen(fileNameSkim,"w");
  FILE *f6 = fopen(fileNameChist,"w");
  FILE *f7 = fopen(fileNameSkimnps,"w");

// open the rootfile to be analyzed
  TString fileNameD = 
    //"/cache/hallc/c-nps/analysis/online/replays/production/" ;
"/cache/hallc/c-nps/analysis/pass1/replays/skim/" ;
  TString fileNameD2nd = 
"/cache/hallc/c-nps/analysis/online/replays/production/" ;
 fileNameD += "nps_hms_skim_"; //read the root file from data
 fileNameD2nd += "nps_hms_skim_"; //read the root file from data
 fileNameD += runNumber; //read the root file from data
 fileNameD2nd += runNumber; //read the root file from data
 fileNameD += "_1_-1.root"; //read the root file from data
 fileNameD2nd += "_1_-1.root_2nd"; //read the root file from data
 
  Int_t foundfile =0 ;
  TFile *f1 = new TFile(fileNameD);
  if(f1->GetSize()!=-1) foundfile = 2 ; 
  if(foundfile == 2) {
    
  TTree *tt = (TTree*)f1->Get("T");
  //get the relevant branch
  int nentriesD = tt->GetEntries();
  cout<<"Entries:\t"<<nentriesD<<endl;

  gROOT->SetBatch(kTRUE);
  
  Double_t HgtrX, HgtrTh, HgtrPh, hdelta ; 
  Double_t hdcx,hdcxp,hdcy,hdcyp ; 
  Double_t evtType, HgtrP ;
  Double_t cointime, HhodStatus, starttime, fptime, paeronpe,paeropt[99],paeront[99],paeroptd[99],paerontd[99];
  Double_t ntime,xexit,yexit ; 
  Double_t hbeta, hcalepr, hcaletot, hcernpe ;
  Double_t sum1,sum2,counters[10],ethcntr[10];
  Double_t  helmpsa,helnega,helposa,hgtry,pgtry,pcaltot,pgdsc,gdtrk,trkeff,trkeffer,prf, hrf, prfraw, hrfraw, hztar, pztar,trkeffg,trkeffger,trkefft,trkeffter,htof1,htof2 ;
  Double_t  hgdsc,htrkeff,htrkeffer,phbig,ppbig,t1big,t2big,sum3,sum4,sum5,sum6,sum7,sum8 ;
  Double_t hel1,helpred,helrep,helmps,helnqrt,hcaltot ;  
  Int_t  helmpsi,helnegi,helposi,icc,chist[100],chistp[100],chistpt[100];
  Int_t chistpg[100],chistt[100],chistg[100], cltimehist6[100] ;
  Int_t t16diffh[100],td1h[100],td6h[100], td16h[100] ;
  Int_t ctxh[100][20], ctyh[100][20],eclhist[100][20];
  Int_t cteh[100][20] ;
  Int_t ctimehist[120], cltimehist[120], vcltimehist[100] ;
  Int_t evtypecnr[10],ptdcmulth[15],icluster;
  Double_t trig1tdc,trig6tdc,edtmtdc ;
  Double_t  hdchit1,hdchit2,hntrk,ngcorr ;
  Double_t  clusE[10000],clusX[10000],clusY[10000],clusT[10000] ;
  Double_t  vclusE[10000],vclusX[10000],vclusY[10000],vclusT[10000] ;
  Int_t vclusSize[10000],goodcl[10000] ;
  int vnclus ;
  Double_t nclust,block_e[2000],cluster_ID[2000] ;
  Double_t block_x[2000],block_y[2000],block_t[2000] ;
  Double_t ctpi1,ctpi2,hst,pst,hpst,avm[30];
  Double_t avmb[1080],sumb[1080],minion[1080],minioner[1080];
  Int_t nchist[100], mallhist[100], mbesthist[100] ;
  Int_t cerhist[100], ephist[100],nc1hist[100],nc2hist[100] ; 
  Int_t epihist[100], mvbesthist[100], msbesthist[100] ;
  Int_t cthist1[120],cthist6[120],cthist16[120],cthist0[120] ;  
  Int_t eoverphist[100],blktimeh[1080][20] ;  
  Int_t blkmassh[1080][20],mbycolhist[100][30] ;  
  Int_t eloblkhist[1080][20], ehiblkhist[1080][20] ;
  Int_t nevnt=0, trg1,trg6,ic1,ic2,iclo,ichi,clth[100][5];
  Int_t nevntall=0;
  Double_t enorm = 30. ; 
  Double_t emin = 0.7 ; 

  if(runNumber>1632) enorm = 2.4 ;
  if(runNumber>1680) enorm = 1.35 ;
// this is slope for enorm versus column
//    double eslope = -0.13 / 60. ;
// use this to calibrate individual blocks
  if(runNumber>0) enorm = 1.35 * 0.93 ;
// new calibration with ep elastic starting here
  if(runNumber>5236) enorm = 1.095 ;
// for pass1, hope this works:
  enorm = 1.000 ;
    double eslope = 0. ;

  for(icc=0 ; icc<100 ;  icc++) {
    nchist[icc]=0. ;
    nc1hist[icc]=0. ;
    nc2hist[icc]=0. ;
    cerhist[icc]=0. ;
    ephist[icc]=0. ;
    mallhist[icc]=0. ;
    mbesthist[icc]=0. ;
    mvbesthist[icc]=0. ;
    msbesthist[icc]=0. ;
    epihist[icc]=0. ;
    ctimehist[icc]=0 ;
    cthist0[icc]=0 ;
    cthist1[icc]=0 ;
    cthist6[icc]=0 ;
    cthist16[icc]=0 ;
    cltimehist[icc]=0 ;
    cltimehist6[icc]=0 ;
    vcltimehist[icc]=0 ;
    t16diffh[icc]=0 ;
    td1h[icc]=0 ;
    td16h[icc]=0 ;
    td6h[icc]=0 ;
    for(int jj=0 ; jj<20 ; jj++) {
      ctxh[icc][jj]=0 ;
      ctyh[icc][jj]=0 ;
      eclhist[icc][jj]=0 ;
      cteh[icc][jj]=0 ;
    }
    for(int jj=0 ; jj<20 ; jj++) {
      mbycolhist[icc][jj]=0 ;
    }
  }
  for(icc=0 ; icc<1080 ;  icc++) {
    for(int jj=0 ; jj<20 ; jj++) {
      blktimeh[icc][jj]=0 ;
      blkmassh[icc][jj]=0 ;
      eloblkhist[icc][jj]=0 ;
      ehiblkhist[icc][jj]=0 ;
    }
  }
  for(icc=0 ; icc<100 ;  icc++) {
    for(int jj=0 ; jj<5 ; jj++) {
      clth[icc][jj]=0 ;
    }
  }
  for(icc=0 ; icc<120 ;  icc++) {
    ctimehist[icc]=0 ;
    cthist0[icc]=0 ;
    cthist1[icc]=0 ;
    cthist6[icc]=0 ;
    cthist16[icc]=0 ;
  }
  for(icc=0 ; icc<10 ;  icc++) {
    counters[icc]=0 ;
    ethcntr[icc]=0 ;
    evtypecnr[icc]=0 ;
  }
 
  tt->SetBranchAddress("H.dc.x_fp",  &hdcx); 
  tt->SetBranchAddress("H.dc.y_fp",  &hdcy); 
  tt->SetBranchAddress("H.dc.xp_fp", &hdcxp); 
  tt->SetBranchAddress("H.dc.yp_fp", &hdcyp); 
  tt->SetBranchAddress("H.gtr.y", &hgtry); 
  tt->SetBranchAddress("H.gtr.x", &HgtrX); 
  tt->SetBranchAddress("H.gtr.p", &HgtrP); 
  tt->SetBranchAddress("H.gtr.beta", &hbeta); 
  tt->SetBranchAddress("H.gtr.dp", &hdelta);   
  tt->SetBranchAddress("H.gtr.th", &HgtrTh);   
  tt->SetBranchAddress("H.gtr.ph", &HgtrPh);   
  tt->SetBranchAddress("H.react.z", &hztar);  
  tt->SetBranchAddress("H.cal.eprtracknorm", &hcalepr);
  tt->SetBranchAddress("H.cal.etottracknorm", &hcaletot);  
  tt->SetBranchAddress("H.cal.etotnorm", &hcaltot); 
  tt->SetBranchAddress("H.cer.npeSum", &hcernpe); 
  tt->SetBranchAddress("H.hod.goodscinhit", &hgdsc);

  tt->SetBranchAddress("T.helicity.hel",&hel1) ;
  tt->SetBranchAddress("T.helicity.helpred",&helpred) ;
  tt->SetBranchAddress("T.helicity.helrep",&helrep) ; 
  tt->SetBranchAddress("T.helicity.mps",&helmps) ;
  tt->SetBranchAddress("T.helicity.nqrt",&helnqrt) ;
                    
  tt->SetBranchAddress("H.hod.starttime", &starttime);       

  tt->SetBranchAddress("CTime.epCoinTime1_ROC1", &ctpi1);
  tt->SetBranchAddress("CTime.epCoinTime2_ROC1", &ctpi2);
  tt->SetBranchAddress("H.hod.fpHitsTime", &fptime); 
  tt->SetBranchAddress("H.dc.ntrack", &hntrk); 
  
  tt->SetBranchAddress("NPS.cal.nclust", &nclust); 
  tt->SetBranchAddress("NPS.cal.clusE", &clusE); 
  tt->SetBranchAddress("NPS.cal.clusX", &clusX); 
  tt->SetBranchAddress("NPS.cal.clusY", &clusY); 
  tt->SetBranchAddress("NPS.cal.clusT", &clusT); 
  tt->SetBranchAddress("NPS.cal.fly.block_clusterID",&cluster_ID);
  tt->SetBranchAddress("NPS.cal.fly.e",&block_e) ;
  tt->SetBranchAddress("NPS.cal.fly.x",&block_x) ;
  tt->SetBranchAddress("NPS.cal.fly.y",&block_y) ;
  tt->SetBranchAddress("NPS.cal.fly.goodAdcTdcDiffTime",&block_t) ;
  tt->SetBranchAddress("NPS.cal.vtpClusE",&vclusE);
  tt->SetBranchAddress("NPS.cal.vtpClusSize",&vclusSize);
  tt->SetBranchAddress("NPS.cal.vtpClusTime",&vclusT);
  tt->SetBranchAddress("Ndata.NPS.cal.vtpClusY",&vnclus);
  tt->SetBranchAddress("T.hms.npsTRIG1_tdcTimeRaw",&trig1tdc);
  tt->SetBranchAddress("T.hms.npsTRIG6_tdcTimeRaw",&trig6tdc);
  tt->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw",&edtmtdc);
  Double_t ctimepi, ctimepinew, ctimeK, ctimep,ctimeraw,tmp ;                          
//   fprintf(f2,"  0.00  0.00 0.00\n") ;

// Start of loop over events
  for (int kk=0; kk<nentriesD;  kk++){
    //    for (int kk=0; kk<60000;  kk++){
    if((kk/100000)*100000 == kk) printf("doing kk=%d run=%d\n",kk,runNumber) ;
    tt->GetEntry(kk);
    evtType = tt->GetLeaf("fEvtHdr.fEvtType")->GetValue(); 
    if(kk>0 && kk<3) printf("kk,evtype %d %3.1f edtm=%8.1f %8.1f %8.1f \n",kk,evtType,edtmtdc,trig1tdc, trig6tdc) ;

// normalize cluster energies. Also get highest energy cluster
    int nhi = -1 ; 
    double ehi = 0. ;
    for(int n1=0 ; n1 < nclust ; n1++) {
      clusE[n1] *= enorm * (1. + eslope * (clusX[n1] + 30.)) ;
      if(clusE[n1] > ehi) {
	nhi = n1 ;
	ehi = clusE[n1] ;
      }
    }

    Double_t timeoffset = -0.1 ;
    Int_t irun = runNumber ;
    if(irun >= 3057 && irun < 3700 ) timeoffset = 4.1 ;
    if(irun >= 3700 && irun < 4634 ) timeoffset = 0.8 ;
    if(irun >= 4634 && irun < 4965 ) timeoffset = -5.5 ;
    if(irun >= 4965 && irun < 5350 ) timeoffset = 1.0 ;
    if(irun >= 5350 && irun < 5523 ) timeoffset = 4.4 ;
    if(irun >= 5523 && irun < 9999 ) timeoffset = 2.9 ;

    int nclusgood = 0 ;
    int nintime = 0 ;
    if(nclust > 0){
      for(int i=0 ; i < 1080 ; i++) {
	block_t[i] = block_t[i]  + timeoffset ; ;
      }
      for(int n1=0 ; n1 < nclust ; n1++) {
	goodcl[n1]=0 ;
        if(clusE[n1] > -999. && clusE[n1] <999. &&
           clusX[n1] > -999. && clusX[n1] <999. &&
           clusY[n1] > -999. && clusY[n1] <999. &&
           clusT[n1] > 0 && clusT[n1 ] < 200.) {
          clusT[n1] = clusT[n1] + timeoffset ;
	  int it = int(clusT[n1] - 100.) ;
	  cthist0[it]++ ;
	  if(clusE[n1] > 0.6) cthist1[it]++ ;
	  if(clusE[n1] > 1.0) cthist6[it]++ ;
	  if(clusE[n1] > 0.6 && clusT[n1]>130. &&
	     clusT[n1] < 170.) {
	    nclusgood++ ;
	    goodcl[n1]=1 ;
	  }
	  if(clusE[n1] > 0.6 && clusT[n1]>145. &&
	     clusT[n1] < 155.) nintime++ ;
	}
      }
    }

// ncluster histograms
    icc = nclust ; 
    if(icc < 0) icc = 0 ;
    if(icc > 99) icc = 99 ;
    nchist[icc]++ ;
    icc = nclusgood ; 
    if(icc < 0) icc = 0 ;
    if(icc > 99) icc = 99 ;
    nc2hist[icc]++ ;
    icc = nintime ; 
    if(icc < 0) icc = 0 ;
    if(icc > 99) icc = 99 ;
    nc1hist[icc]++ ;


    counters[0] = counters[0]+1 ;

    if(edtmtdc<0.1 && hdelta > -12. && hdelta < 12.) {
      counters[1] = counters[1]+1 ;
    }
    if(edtmtdc<0.1 && hdelta > -12. && hdelta < 12. &&
       hcaltot > 0.6) {
      counters[2] = counters[2]+1 ;
    }
    if(edtmtdc<0.1 && hdelta > -12. && hdelta < 12. &&
       hcaltot > 0.6 && hcernpe > 1.0 ) {
      counters[3] = counters[3]+1 ;
    }
    if(edtmtdc<0.1 && hdelta > -12. && hdelta < 12. &&
       hcaltot > 0.6 && hcernpe > 1.0 && nclusgood>1) {
      counters[4] = counters[4]+1 ;
    }
    if(edtmtdc<0.1 && hdelta > -12. && hdelta < 12. &&
       hcaltot > 0.6 && hcernpe > 1.0 && nclusgood>1
       && nintime>0) {
      counters[5] = counters[5]+1 ;
    }
// HMS  Cerenkov
    if(edtmtdc<0.1 && hdelta > -12. && hdelta < 12. &&
       hcaltot > 0.6 && nclusgood>1
       && nintime>0 && ehi>1.0) {
      if(hcernpe >= 0 && hcernpe <20) {
	int j = (hcernpe * 5) ;
	cerhist[j]++ ;
      }
    }
// HMS  E/p
    if(edtmtdc<0.1 && hdelta > -12. && hdelta < 12. &&
       hcernpe > 1.0 && nclusgood>1
       && nintime>0 && ehi>1.0) {
      if(hcaltot >= 0 && hcaltot <2.) {
	int j = (hcaltot * 50) ;
	ephist[j]++ ;
      }
    }
// skip etdm events and no-track HMS events and also
// require good e/p snf cer (FIX if doing protons!)
// also require two cluster E>0.6
// and at least one cluster 145<t<155
// add highest energy cluster over 1 gev
    if(edtmtdc<0.1 && hdelta > -12. && hdelta < 12. &&
       hcaltot > 0.6 && hcernpe > 1.0 && nclusgood>1
       && nintime>0 && ehi>1.0) {
      counters[6] = counters[6]+1 ;
 
// Get helicity: differene for sping18 and later
// xxx need to fix
      helmpsi = 0 ;
      if(helmps > 0) helmpsi = 1 ;
      helnegi = 0 ;
      if(hel1 < 0) helnegi = 1 ;
      helposi = 0 ;
      if(hel1 > 0) helposi = 1 ;

// coin time histogram all events
      icc = ctpi1 - 50 ;
      if(icc < 0) icc = 0 ;
      if(icc > 119) icc = 119 ;
      ctimehist[icc]++ ;

// cluster time hist all events
      for(int n1=0 ; n1 < nclust ; n1++) {
	icc = (int)clusT[n1] - 105 ; 
	if(icc < 0) icc = 0 ;
	if(icc > 99) icc = 99 ;
	if(clusE[n1]>0.3) cltimehist[icc]++ ;
      }

      counters[2] = counters[2]+1 ;

      trg1 = 0 ;
      if(trig1tdc>31300 && trig1tdc<35000) trg1 = 1;
      trg6 = 0 ;
      if(trig6tdc>31300 && trig6tdc<35000) trg6 = 1;

// M2g mass spectrum
      Double_t mbest = 0. ;
      if(nclusgood == 2) {
	Double_t dnps = 307. ; // targ to NPS dist.
	if(runNumber > 3700) dnps = 407. ; // 4 m in 2024
	if(runNumber > 4633) dnps = 607. ; // 6 m in 2024 for small angle NPS
	if(runNumber > 4965) dnps = 407. ; // back to 4 m on 3/4/2024
	if(runNumber > 5350) dnps = 307. ; // back to 3 m on 3/19/2024
	if(runNumber > 5522) dnps = 357. ; // went to 3.5 m 
	for(int n1=0 ; n1 < nclust-1 ; n1++) {
	  if(goodcl[n1]) {
	    for(int n2=n1+1 ; n2 < nclust ; n2++) {
	      if(goodcl[n2]) {
		Double_t tdiff = clusT[n1] - clusT[n2] ;
		if(tdiff > -3. && tdiff < 3.) {
		  Double_t xdiff = clusX[n1] - clusX[n2] ;
		  Double_t ydiff = clusY[n1] - clusY[n2] ;
		  Double_t theta = sqrt(xdiff * xdiff + ydiff * ydiff) / dnps ;
		  Double_t sth2 = sin(theta/2.) ;
		  Double_t sinsq = sth2 * sth2 ; 
		  Double_t mass = sqrt(4. * clusE[n1] * clusE[n2] * sinsq) ;
		  Double_t chkk = sqrt(2.*clusE[n1]*clusE[n2]*theta*theta/2.) ;
		  mbest = mass ;
		}
	      }
	    }
	  }
	}
	icc = (int)(mbest / 0.135 * 50.) ;
	if(icc < 0) icc = 0 ;
	if(icc > 99) icc = 99 ;
	mbesthist[icc]++ ;
      }
// writing to skim file here
      fprintf(f2,
"%3d %1d %1d %1d %6.2f %7.4f %7.4f %6.2f %7.4f %7.4f %6.2f %6.2f %7.4f %7.4f %7.3f   \n",
	  nclusgood,(int)hel1+1,(int)helrep+1,(int)helmps+1,
          hdelta,HgtrTh,HgtrPh,hcernpe,hcalepr,hcaletot,
	  hdcx, hdcy, hdcxp, hdcyp,hztar) ;
      for(int n1=0 ; n1 < nclust ; n1++) {
	if(goodcl[n1]>0) {
	  int nn = 0 ;
	  for(int i=0 ; i < 1080 ; i++) {
	    icluster = (int)cluster_ID[i] ; 
	    if(icluster > n1 - 0.5  && icluster < n1+0.5
	       && block_e[i]>0.) nn++ ;
	  }
	  fprintf(f2,"%4d %7.2f %7.2f %7.2f %7.2f\n",
		  nn,clusE[n1],clusX[n1],clusY[n1],clusT[n1]) ;
	  for(int i=0 ; i < 1080 ; i++) {
	    icluster = (int)cluster_ID[i] ; 
	    if(icluster > n1 - 0.5  && icluster < n1+0.5
	       && block_e[i]>0.) {
	      fprintf(f2,"%4d %3d %8.4f %8.1f \n",
		      i,icluster,block_e[i],block_t[i]) ;
              if(block_e[i] > 0.3 * clusE[n1]) {
		int jj = block_t[i] - 140. ;
		if(jj > -1 && jj<20) blktimeh[i][jj]++ ;
		jj = (mbest - 0.115) * 20. / 0.040 ;
		if(jj > -1 && jj<20) blkmassh[i][jj]++ ;
	      }
	    }
	  }
	}
      }
    }  // edtm check
  } // loop over events
  
// write out histograms
  for(icc=0 ; icc<100 ; icc++){
    printf("%2d %7d %6d %6d %6d %6d %6d %6d %6d %4d \n",
      icc,nchist[icc],nc2hist[icc],nc1hist[icc],
      cerhist[icc],ephist[icc],
      cthist0[icc],cthist1[icc],cthist6[icc],
	   mbesthist[icc]);
    fprintf(f6,"%2d %7d %6d %6d %6d %6d %6d %6d %6d %4d \n",
      icc,nchist[icc],nc2hist[icc],nc1hist[icc],
      cerhist[icc],ephist[icc],
      cthist0[icc],cthist1[icc],cthist6[icc],
	   mbesthist[icc]);
  }

  for(icc=0 ; icc<8 ; icc++){
    printf("%d  %7.0f %7.3f \n",icc,counters[icc],
	   counters[icc]/counters[0]) ;
    fprintf(f6,"%d  %7.0f %7.3f \n",icc,counters[icc],
	   counters[icc]/counters[0]) ;
  }


// by block outputs
  for(int i=0 ; i<1080 ; i++){
    fprintf(f7,"%4d ",i) ;
    for(int j=0 ; j<19 ; j++) fprintf(f7," %2d",blktimeh[i][j]) ;
    fprintf(f7," %2d \n",blktimeh[i][19]) ;
  }
  for(int i=0 ; i<1080 ; i++){
    avmb[i]=0. ;
    sumb[i]=0. ;
    for(int j=0 ; j<20 ; j++) {
      double mm = 0.115 + 0.040 * (j+0.5)/20. ;
      avmb[i] = avmb[i] + mm * blkmassh[i][j] ; 
      sumb[i] += blkmassh[i][j] ; 
    }
    fprintf(f7,"%4d ",i) ;
    for(int j=0 ; j<19 ; j++) 
      fprintf(f7," %2d",blkmassh[i][j]) ;
    fprintf(f7," %2d \n",blkmassh[i][19]) ;
  }


  fprintf(f2,"  0.00  0.00 0.00\n") ;

  fclose(f2) ;
  fclose(f6) ;
  fclose(f7) ;
  } else {
      printf("--------\n") ;
      printf("run %d not FOUND\n",runNumber) ;
      printf("--------\n") ;

    gSystem->Exit(1) ;
  } 

  gROOT
->SetBatch(kFALSE);

  gSystem->Exit(1) ; 
}

