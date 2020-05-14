#ifndef my_fit_h
#define my_fit_h

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TFile.h"
#include "TCanvas.h"



using namespace std;

float** MyFIT_MGauss(TH1F *data, double *init_params, int ng,float lowb, float upb, bool draw_opt);
TH1F* ReadData(string filename, string nADC);
TH1F* ReadFileList(string filelist, string nADC);

#endif
