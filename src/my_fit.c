#include "my_fit.h"
#include "functions.h"

float** MyFIT_MGauss(TH1F *data, double *init_params, int ng,float lowb, float upb, bool draw_opt){

  /*
    Function to Fit Sum of Gaussian Curves (1 - 5), 
    ng = number of gauss
    lowb, upb = lower and upper bounds of fit
    draw_opt  = Plot on/off (to draw result need a previously opened canvas)
   */

  int n_par = ng * 3;   //Number of Parameters
  
  TF1 *fitFcn;          

  //............ Setting the Correct fit function .....................
  
  if(ng == 1) fitFcn = new TF1("fitFcn",Gauss, lowb, upb, n_par);
  if(ng == 2) fitFcn = new TF1("fitFcn",GaussSum2, lowb, upb, n_par);
  if(ng == 3) fitFcn = new TF1("fitFcn",GaussSum3, lowb, upb, n_par);
  if(ng == 4) fitFcn = new TF1("fitFcn",GaussSum4, lowb, upb, n_par);
  if(ng == 5) fitFcn = new TF1("fitFcn",GaussSum5, lowb, upb, n_par);

  //...........Setting initial Parameters..............................
  
  fitFcn->SetParameters(init_params);

  //.................. Starting The Fit Procedure ...................
  
  if(ng == 1)
    cout<<"\n----------------------   Single Gaussian Fit  ----------------------------\n";
  else if (ng > 1)
    cout<<"\n----------------------  Multiple Gaussian Fit ----------------------------\n";
  
  
  data->Fit("fitFcn","0","",lowb, upb);

  // ....  Calulating ChiSquare  ....
  
  TF1 *fitres = data->GetFunction("fitFcn");
  double Chi = fitres->GetChisquare();
  double NDF = fitres->GetNDF();
  cout<<"CHI QUADRO = "<<Chi/NDF<<endl;

  Double_t par[n_par];
  const Double_t * parErr;
  fitFcn->GetParameters(par);
  parErr = fitFcn->GetParErrors();
 
  //.................. Drawing Results (Requires a Canvas) .....................................
  if(draw_opt){
    
    
      //plot data
      data->Draw();

      fitFcn->SetNpx(5000);
      fitFcn->SetLineWidth(2);
      fitFcn->SetLineColor(kRed);

      //build up legend
      TLegend *legend=new TLegend(0.75,0.75,0.99,0.99);
      legend->SetTextFont(32);
      legend->SetTextSize(0.04);
      legend->AddEntry(data,"Data","l");
      legend->AddEntry(fitFcn,"Fit Result","l");
    
      TF1 *gauss[ng];
    
      for(int i = 0; i < ng; i ++){
      
	string gname = "gauss_"+to_string(i);
	gauss[i] = new TF1(gname.c_str(),Gauss,0,10000,3);
	gauss[i]->SetParameters(&par[i*3]);
	gauss[i]->SetLineColor(kBlue);
	gauss[i]->SetNpx(500000);
	gauss[i]->Draw("same");

	if(ng == 1){
	  legend->AddEntry((TObject*)0,Form("mean = %.2f",par[0]), "");
	  legend->AddEntry((TObject*)0,Form("sigma = %.2f", par[1]), "");
	}else
	  legend->AddEntry((TObject*)0,Form("mean _%d = %.2f",i,par[i*3]), "");
	  legend->AddEntry((TObject*)0,Form("sigma_%d = %.2f", i, par[(1+(i*3))]), "");
      
      }
    
   
      fitFcn->Draw("same");
 
      legend->AddEntry(gauss[0],"Full Gauss","l");
      legend->Draw();
    
  }

  float ** out_par = new float*[n_par];
  for(int i=0 ; i<n_par ; i++) out_par[i] = new float[2];
  for(int i=0 ; i<n_par ; i++){
    out_par[i][0] = par[i];
    out_par[i][1] = parErr[i];
  }
  
  return out_par;
};





TH1F* ReadData(string filename, string nADC){

  filename = "Data/"+filename;
  TH1F *data,*gdata;
  
  
  TFile *fin   = new TFile(Form("%s", filename.c_str()), "read");
  cout <<"    File "<< filename <<" Opened.... "<<endl;
  
  string cname = "cDumpEnergy_"+nADC;  
  TCanvas *c   = (TCanvas*)  fin->Get(cname.c_str());                 
  data   = (TH1F*)c->GetListOfPrimitives()->FindObject("hbad");
  gdata   = (TH1F*)c->GetListOfPrimitives()->FindObject("hgood");
  data->Add(gdata);
  cout <<"    Collected Data.... "<<endl;
  
  fin->Close();

 

  data->SetLineColor(kBlack);
  
  return data;

};


TH1F* ReadFileList(string filelist, string nADC){

  TH1F * data, * tmpdata;
  ifstream File(filelist.c_str());
  string tmpstr;
  int nf=0;
  
  if(File)                                                                                                                    
    {                                                                                                                         
      cout << "    ReadFileList - " << filelist << " opened"<< endl;                                                
      while(getline(File, tmpstr)){      //Leggo i valori separati dalle virgole
	cout << "    ReadFileList - Reading " << tmpstr << endl;
	
	if(nf==0){
	  data=ReadData(tmpstr,nADC);
	}else{
	  tmpdata=ReadData(tmpstr,nADC);	  
	  data->Add(tmpdata);
	}
	cout<<"                        ....Done!"<<endl;

	nf++;
      } 

      cout << "\n    ReadFileList - " << nf << " files merged. \n"<< endl;                                                
    }                                                                                                                         
  else{  // if the file doesn't exist                                                                                         
    cerr << "    ReadFileList - " << filelist << " not found. Exit." << std::endl;                                            
    exit(1);                                                                                                                  
  }                                                                                                                           
                                                                                                                              
  File.close();                                                      
  return data;                                                                                                                 
}        

