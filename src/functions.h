#ifndef functions_h
#define functions_h



#include "Rtypes.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

float * Config(string cfgname);

Double_t Gauss(Double_t *x, Double_t *par);

Double_t GaussSum2(Double_t *x, Double_t *par);

Double_t GaussSum3(Double_t *x, Double_t *par);

Double_t GaussSum4(Double_t *x, Double_t *par);

Double_t GaussSum5(Double_t *x, Double_t *par);

Double_t background(Double_t *x, Double_t *par);
  
Double_t GaussSum_Const(Double_t *x, Double_t *par);

float * PeakDistance(float p1[2], float p2[2], float s[2]);

float * FigureOfMerit(float m1[2], float m2[2], float s1[2], float s2[2]);

float GaussFractionPostThr(float thr, int side, double * par);

#endif
