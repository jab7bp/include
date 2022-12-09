#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <math.h>
#include <algorithm>
#include <string>
#include <chrono>

using namespace std::chrono;
#include "/w/halla-scshelf2102/sbs/jboyd/include/include_files.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/GEM_lookups.h"
#include "/w/halla-scshelf2102/sbs/jboyd/include/beam_variables.h"

Double_t ellipse_fcn(Double_t x, Double_t y,
                     Double_t x0, Double_t y0,
                     Double_t a, Double_t b,
                     Double_t theta) // (in degrees)
	{
		double v = 9999.9;
		if( (a == 0.0) || (b == 0.0) ) return v; // Just a precaution
		//Shift the center:
		x -= x0;
		y -= y0;

		// Un-Rotate the axes:
		theta *= TMath::Pi()/180.0; // degrees->radians
		v = x;
		x = x*std::cos(theta) + y*std::sin(theta);
		y = y*std::cos(theta) - v*std::sin(theta);

		// "Scale" the axes
		x /= a;
		y /= b;

		// Calculate the "normalized distance"
		v = x*x + y*y;
		v = std::sqrt(v);

		return v;
	}

Double_t ellipse_fcn(const Double_t *x, const Double_t *par2D)
{
  return ellipse_fcn(x[0], x[1], // "x", "y"
                     par2D[0], par2D[1], // "x0", "y0"
                     par2D[2], par2D[3], // "a", "b"
                     par2D[4]); // "theta" (in degrees)
}

//Run info and lookups
int runnum = 13479;

void test_2D_fit(){
	TFile *TF = new TFile(Form("%i_dxdy.root", runnum), "READ");
  	hdxdy = static_cast<TH2D*>(TF->Get("h_dxdy"));

  	TF2 *ell_fit = new TF2("ell_fit", ellipse_fcn, 0, 3, -0.5, 0.5, 5);
  	hdxdy->Fit("ell_fit");

  	TCanvas *dxdy = new TCanvas("dxdy", "dxdy", 600, 500);
  	hdxdy->Draw("colz");
  	ell_fit->Draw("same");

}
