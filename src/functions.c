#include "functions.h"


//////////////////////////        LETTORE DI FILE .CFG        ///////////////////////////////////////////////

float * Config(string cfgname){

  vector<float> CFG;
  vector<string> TYPES;
  string str, CFGNAME;
  int ncfg=0;

  CFGNAME = "cfg/"+cfgname+".cfg";

  ifstream file(CFGNAME); //Apro il file
  
  // cout<<"leggo... "<<wn<<"\n";
  
  if(!getline(file, str))  cout<<"PROBLEMA DI LETTURA"<<endl;
  
   
  while(getline(file, str)){      //Leggo i valori separati dalle virgole
    string type, val;
    type = str.substr(0,str.find(':')); 
    val = str.substr(str.find(':')+1, string::npos);
    //cout<<"tempo: "<<tim<<"   ampiezza:"<<amp<<endl;
    TYPES.push_back(type);
    CFG.push_back(stof(val));

    cout<<TYPES.at(ncfg)<<":  "<<CFG.at(ncfg)<<endl;
    
    ncfg++;
  }

  //cout<<"ho letto "<<nsamp<<"righe "<<"\n";
  
  float * cfg;
  cfg = new float[ncfg];
  
  for(int i=0; i<ncfg; i++) cfg[i] = CFG.at(i);

  cout<<"n configs: "<<ncfg<<endl;

  return cfg;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////          FITTING   FUNCTIONS     ////////////////////////////////////////////

Double_t Gauss(Double_t *x, Double_t *par) {
  return (par[2]/(par[1]*sqrt(2*TMath::Pi()))*(exp(-(pow(x[0]-par[0],2))*(0.5/(pow(par[1],2))))));
}

Double_t GaussSum2(Double_t *x, Double_t *par) {
  return Gauss(x,par) + Gauss(x,&par[3]);
}

Double_t GaussSum3(Double_t *x, Double_t *par) {
  return Gauss(x,par) + Gauss(x,&par[3])+ Gauss(x,&par[6]);
}

Double_t GaussSum4(Double_t *x, Double_t *par) {
  return Gauss(x,par) + Gauss(x,&par[3])+ Gauss(x,&par[6])+ Gauss(x,&par[9]);
}

Double_t GaussSum5(Double_t *x, Double_t *par) {
  return Gauss(x,par) + Gauss(x,&par[3])+ Gauss(x,&par[6]) + Gauss(x,&par[9]) + Gauss(x,&par[12]);
}

Double_t background(Double_t *x, Double_t *par) {
   return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}


Double_t GaussSum_Const(Double_t *x, Double_t *par){
  return Gauss(x,par) + Gauss(x,&par[3])+ Gauss(x,&par[6]) + Gauss(x,&par[9]) + background(x,&par[12]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////


float * PeakDistance(float p1[2], float p2[2], float s[2]){
  
  float * res = new float[2];
  
  res[0]  = (p1[0]-p2[0])/s[0];
  
  res[1]  = 0;
  res[1] += (p1[1]*p1[1])/(s[0]*s[0]);
  res[1] += (p2[1]*p2[1])/(s[0]*s[0]);
  res[1] += (p1[0]-p2[0])*(p1[0]-p2[0])/((s[0]*s[0])*(s[0]*s[0]))*(s[1]*s[1]);
  res[1]  = sqrt(res[1]);
    
  return res;
};




float * FigureOfMerit(float m1[2], float m2[2], float s1[2], float s2[2]){

  float * res = new float[2];

  float num = abs(m1[0]-m2[0]), den = sqrt(s1[0]*s1[0] + s2[0]*s2[0]); 

  res[0] = num/den;

  //correggere propagazione errore
  res[1]  = 0;
  res[1] += (m1[1]*m1[1])/(den*den);
  res[1] += (m2[1]*m2[1])/(den*den);
  res[1] += (num/(den*den))*(num/(den*den))*(s1[1]*s1[1]);
  res[1] += (num/(den*den))*(num/(den*den))*(s2[1]*s2[1]);
  res[1]  = sqrt(res[1]);

  return res;
};




float GaussFractionPostThr(float thr, int side, double * par){
  /*
    side=0 -> left of thr,  side=1 -> right of thr
  */
  

  int npoints=1000000;
  double* points=new double[npoints]();
  double* gauss=new double[npoints]();

  float tail=0, total_sum=0, ratio=0;
  float x;
  float intmin=par[0] - 20*par[1], intmax=par[0] + 20*par[1];

  for(int i=0;i<npoints;i++){
    x = intmin+((intmax-intmin)/npoints)*i;  
    points[i]=x;
    gauss[i]=Gauss(&points[i], par);
    //  fprintf(stderr,"xvalue = %.3f ,  gauss = %.2f \n",points[i],gauss[i]); 
    total_sum +=gauss[i];
    
    if(side == 1)  if(x>=thr) tail+=gauss[i];
    if(side == 0)  if(x<=thr) tail+=gauss[i];
    //   fprintf(stderr,"x = %.3f ,  tot sum= %.3f,  tail = %.3f\n",x,total_sum,tail);
  }
  
  ratio=tail/total_sum*100;

  delete[] points;
  delete[] gauss;
  
  return ratio;
};
