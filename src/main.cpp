#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>



#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TTree.h> 
#include <TF1.h>
#include <TH1.h>
#include <TLegend.h>
#include <TPaveText.h>
#include "src/functions.h"
#include "src/my_fit.h"




using namespace std;


void reco(string run, string adcn, string cfgfile, int Nfile){
  /*
  //open root file
  string wn    = run + adcn + "-Monitoring"+"RawEnergy.root";
  TFile *fin   = new TFile(Form("%s", wn.c_str()), "read");

  //read data from input file
  string cname = "cDumpEnergy_"+adcn;  
  TCanvas *c   = (TCanvas*)  fin->Get(cname.c_str());                 
  TH1F *data   = (TH1F*)c->GetListOfPrimitives()->FindObject("hbad"); 
  data->SetLineColor(kBlack);
  fin->Close();
  */
  
  TH1F * data,* ndata;
  if(Nfile>1)  data = ReadFileList(run, adcn);
  else   data = ReadData(run, adcn);

  
  //create output file
  string oname = "Spectrum_" + run;
  if(Nfile>1) oname = oname + ".root";
  TFile *fout  = new TFile(Form("Spect/%s", oname.c_str()), "recreate");   
 
  //READ CONFIG FILE
  float *cfg;
  cfg = Config(cfgfile);
  
  //cout<<"ho aperto tutto\n";

  /// TITLE OF THE DATA
  string dname = "Energy Specrum - ADC_"+adcn;
  data->SetTitle(dname.c_str());
  data->Rebin(4);

  //RESCALING TO 1PE
  float spe=cfg[24];     //ch-19: 3040    ch-20: 3840    ch-21: 4040 
  int Nbin, NEntries;
  float maxbin, binrange;
  Nbin=data->GetNbinsX();
  maxbin=data->GetXaxis()->GetXmax();
  NEntries=data->GetEntries();
  binrange=maxbin/spe;
  ndata = new TH1F("ndata",dname.c_str(),Nbin,0,binrange);
  float bin_v;
  for(int i=0; i<Nbin; i++){
    bin_v=data->GetBinContent(i);
    ndata->SetBinContent(i,bin_v);
  }
  ndata->SetLineColor(kBlack);
  //ndata->Scale(NEntries);


  //---------------- Rate Histo
  float norm_r=cfg[25];  //49980.019-maj4    42409.990-maj1   90511.669-maj1_l  653.35328-piedi 
  float cvalue,dr;
  string rname=" Photoelectron Rate (ch-"+adcn+") ; Threshold(PE) ; rate (Hz)";
  TH1F * hrate = new TH1F("hrate",rname.c_str(),Nbin,0,binrange);
  hrate->SetBinContent(0,NEntries/norm_r);
  for(int i=1; i<Nbin;i++){
    dr=data->GetBinContent(i)/norm_r;
    cvalue=hrate->GetBinContent(i-1);
    hrate->SetBinContent(i,cvalue-dr);
  }
  // hrate->Scale(1);
  hrate->SetLineColor(kBlack);
  hrate->Write();
  
  //-----------------------------

  
  //  data->Scale(1/spe);

  //                                                                             //
  //   ____________________________               ____________________________   //
  //   ____________________________    F  I  T    ____________________________   //
  //                                                                             //
  
  //               Initialize TFits, Setting ranges and parameters               //
 
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  int n_par1 = 3*cfg[0];                                   //n parameter
  float rf_min1 = cfg[1], rf_max1 = cfg[2];                //fit range

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  int n_par2 = 3*cfg[3];                                   //n parameter
  float rf_min2 = cfg[4], rf_max2 = cfg[5];                //fit range
  
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  int n_par3 = 3*cfg[6];                                   //n parameter
  float rf_min3 = cfg[7], rf_max3 = cfg[8];                //fit range
  
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  

  //   +-------------------------------------------------------------+
  //   |                                                             |
  //   |     +-----------------------------------------------+       |
  //   |     |     SETTING INITIAL PARAMETERS FOR FIT        |       |
  //   |     +-----------------------------------------------+       |
  //   |                                                             |
  //   +-------------------------------------------------------------+
  
  Double_t* init_params;
  init_params= new Double_t[15];

  //----------------------------------------------------------------+
  init_params[0]  = cfg[9];   //µ1
  init_params[1]  = cfg[10];  //s1           NOISE PEAK
  init_params[2]  = cfg[11];  //norm1
  //----------------------------------------------------------------+
  init_params[3]  = cfg[12];  //µ2
  init_params[4]  = cfg[13];  //s2          FIRST PE PEAK
  init_params[5]  = cfg[14];  //norm2
  //----------------------------------------------------------------+
  init_params[6]  = cfg[15];  //µ3
  init_params[7]  = cfg[16];  //s3          SECOND PE PEAK
  init_params[8]  = cfg[17];  //norm3 
  //----------------------------------------------------------------+
  init_params[9]  = cfg[18];  //µ4
  init_params[10] = cfg[19];  //s4          THIRD PE PEAK
  init_params[11] = cfg[20];  //norm4
  //----------------------------------------------------------------+
  init_params[12] = cfg[21];  //µ5
  init_params[13] = cfg[22];  //s5         BACKGROUND STUFF
  init_params[14] = cfg[23];  //norm5
  //----------------------------------------------------------------+
 
  //
  //                    FITTING AND PRINTING RESULTS
  //

  //PREPARING CANVAS AND PADS.....
  TCanvas *c1 = new TCanvas("c1","My Canvas",200,10,700,900);
  TPad *pad1  = new TPad("pad1","Multiple Gaussian Fit", 0.001, 0.501, 0.999, 0.999);
  TPad *pad2  = new TPad("pad2","Single Gaussian Fit"  , 0.001, 0.001, 0.499, 0.499);
  TPad *pad3  = new TPad("pad3","Single Gaussian Fit"  , 0.501, 0.001, 0.999, 0.499);
  
  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  
  //-----------------------------------  PAD 1 -----------------------------------------
  
  pad1->cd();

  float ** out_params1;
  if(cfg[0])
    {
      if(rf_max1>binrange)
	{
	  fprintf(stderr,"Error: bad fitting range in fit 1, skipping this fit.");
	  cfg[0]=0;
	}
      else out_params1 = MyFIT_MGauss(ndata, &init_params[0], cfg[0], rf_min1, rf_max1, true);
    }
  else ndata->Draw();
  
  //-----------------------------------  PAD 2 -----------------------------------------
  
  pad2->cd();

  float ** out_params2;
  if(cfg[3])
    {
     if(rf_max2>binrange)
	{
	  fprintf(stderr,"Error: bad fitting range in fit 3, skipping this fit.");
	  cfg[3]=0;
	}
      else out_params2 = MyFIT_MGauss(ndata, &init_params[3], cfg[3], rf_min2, rf_max2, true);
    }
  else ndata->Draw();
  
  //-----------------------------------  PAD 3 -----------------------------------------
  
  pad3->cd();

  float ** out_params3;
  if(cfg[6])
    {
      if(rf_max3>binrange)
	{
	  fprintf(stderr,"Error: bad fitting range in fit 3, skipping this fit.");
	  cfg[6]=0;
	}
      else out_params3 = MyFIT_MGauss(ndata, &init_params[0], cfg[6], rf_min3, rf_max3, true);
    }
  else ndata->Draw();
  
  //----------------------------------    INTEG   --------------------------------------

  //------------------ ratio 1 + ratio 2
  float ratio=0,ratio2=0;
  
  double* noisepar=new double[3]();
  double* spepar=new double[3]();
    
  if(cfg[0]==2){

    for(int i=0;i<3;i++){
      noisepar[i]=out_params1[i][0];
      spepar[i]=out_params1[i+3][0];
    }
    float thr=0.3;
    ratio  = GaussFractionPostThr(thr, 1, noisepar);
    ratio2 = GaussFractionPostThr(thr, 1, spepar);
    
    
    
    /*
    int npoints=1000000;
    double* points=new double[npoints]();
    double* points_spe=new double[npoints]();
    double* noiseg=new double[npoints]();
    double* speg=new double[npoints]();
    double* noisepar=new double[3]();
    double* spepar=new double[3]();
    for(int i=0;i<3;i++){
      noisepar[i]=out_params1[i][0];
      spepar[i]=out_params1[i+3][0];
    }
    float intmin_spe=out_params1[3][0] - 20*out_params1[4][0],intmax_spe=out_params1[3][0] + 20*out_params1[4][0];
    float intmin=out_params1[0][0] - 10*out_params1[1][0],intmax=out_params1[0][0] + 10*out_params1[1][0];
    float tail=0,tail2=0, total_sum=0,total_sum2=0;


    for(int i=0;i<npoints;i++){
      points[i]=intmin+((intmax-intmin)/npoints)*i;
      points_spe[i]=intmin_spe+((intmax_spe-intmin_spe)/npoints)*i;  
      noiseg[i]=Gauss(&points[i],noisepar);
      speg[i]=Gauss(&points_spe[i],spepar);
      //  fprintf(stderr,"xvalue = %.3f ,  gauss = %.2f \n",points[i],noiseg[i]);
    };

    float thr=(0.3*out_params1[3][0]);
    //    fprintf(stderr,"\n treshold = %.3f \n\n",thr);
  
    for(int i=0;i<npoints;i++){
      float x = intmin+((intmax-intmin)/npoints)*i;
      float x2 = intmin_spe+((intmax_spe-intmin_spe)/npoints)*i;
      total_sum +=noiseg[i];
      total_sum2+=speg[i];
    
      if(x>=thr) tail+=noiseg[i];
      if(x2<=thr) tail2+=speg[i];
      //   fprintf(stderr,"x = %.3f ,  tot sum= %.3f,  tail = %.3f\n",x,total_sum,tail);
    }

    ratio=tail/total_sum*100;
    ratio2=(1 - (tail2/total_sum2))*100;
    */
  }

  
  
  //----------------------------------    PRINT   --------------------------------------
  
  float *fom1=0, *fom2=0;
  if(cfg[0]==2)         fom1 = FigureOfMerit(out_params1[0],  out_params1[3], out_params1[1], out_params1[4]);
  if(cfg[3]==cfg[6]==1) fom2 = FigureOfMerit(out_params2[0],  out_params3[0], out_params2[1], out_params2[1]);
    
  if(cfg[0]==2)
    fprintf(stderr,"\033[1m\033[92m\n  FIGURA DI MERITO (1)   =  %.4f +/- %.4f \033[0m\n\n",fom1[0],fom1[1]);
  if(cfg[3]==cfg[6]==1)
    fprintf(stderr,"\033[1m\033[92m  FIGURA DI MERITO (2)   =  %.4f +/- %.4f \033[0m\n\n",fom2[0],fom2[1]);
  fprintf(stderr,"\033[1m\033[92m  RUMORE POST SOGLIA     =  %.4f %%  \033[0m\n\n",ratio);
  fprintf(stderr,"\033[1m\033[92m  SPE POST SOGLIA        =  %.4f %%  \033[0m\n\n",ratio2);
  fprintf(stderr,"\n  NUMERO ENTRIES  = %d\n\n",NEntries);
  
  //------------------------------------------------------------------------------------
  
  c1->Write();
  c1->Close();

  //_______________________________________________________________________________
 

  cout<<"chiudo "<<oname<<"\n";
  //data->Write();
  fout->Close();
  
}


int main(int argc, char *argv[]) {
 
  if (argc < 4) {
    cout << "Usage: ./proc  filename   adc_n   cfg   nfiles(optional, default value = 1)" << endl;
    exit(1); 
  }

  string run   = argv[1];
  string adcn  = argv[2];     // wave number
  string cfg   = argv[3];
  int Nfiles   = 1;

  if(argc>4) Nfiles = atoi(argv[4]);
  
  reco(run, adcn, cfg,Nfiles);
  cout<<"uscito"<<endl;
  return 0; 
}
