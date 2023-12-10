#ifndef PLOT_PIMPING_H
#define PLOT_PIMPING_H

void shade_btw_TF1s_w_3rd(TCanvas *c1, TF1 *f1, TF1 *f2, TF1 *f3, double y_Axis_min, double y_Axis_max, TF1 *f1_orig, TF1 *f2_orig, TF1 *f3_orig) {
	//shade the area between f1 and f2 and draw f3 on top

	//create a TGraph to store the function values
	//shaded area is the fill/color/style of f1
	TGraph *gr = new TGraph();
	gr->SetFillColor(f1->GetFillColor());
	gr->SetFillStyle(f1->GetFillStyle());
	f3->Draw("lsame");
	c1->Update();
	//get picture range
	Double_t xmin = f3->GetXmin();
	Double_t xmax = f3->GetXmax();
	Double_t ymin = y_Axis_min;
	Double_t ymax = y_Axis_max;

	//process first function
	Int_t npx = f3->GetNpx();
	Int_t npoints=0;
	Double_t dx = (xmax-xmin)/npx;
	Double_t x = xmin + 0.5*dx;
	while (x <= xmax) {
		Double_t y = f1->Eval(x);
		if (y < ymin) y = ymin;
		if (y > ymax) y = ymax;
		gr->SetPoint(npoints,x,y);
		npoints++;
		x += dx;
	}
	//process second function
	x = xmax - 0.5*dx;
	while (x >= xmin) {
		Double_t y = f2->Eval(x);
		if (y < ymin) y = ymin;
		if (y > ymax) y = ymax;
		gr->SetPoint(npoints,x,y);
		npoints++;
		x -= dx;
	}
	gr->Draw("fsame"); //draw graph with fill area option
	f3->Draw("lsame"); //superimpose function
	f1_orig->Draw("same");
	f2_orig->Draw("same");
}


void fshade_btw_TF1s_w_3rd(TCanvas *c1_orig, TH1D *h_orig_base, TF1 *f1_orig, TF1 *f2_orig, TF1 *f3_orig, double yAxis_min, double yAxis_max) {
	TString h_name = h_orig_base->GetName();
	TString title = h_orig_base->GetTitle();
	TString x_axis_title = h_orig_base->GetXaxis()->GetTitle();
	TString y_axis_title = h_orig_base->GetYaxis()->GetTitle();

	TString canvas_name = Form("%s", h_name.Data());
	canvas_name.ReplaceAll("h_", "c_");

	TString canvas_title = Form("%s", title.Data());

	double xMin = f3_orig->GetXmin();
	double xMax = f3_orig->GetXmax();

	TCanvas *c1 = new TCanvas(canvas_name.Data(), canvas_name.Data(), 600, 500);
	c1->Range(xMin, yAxis_min, xMax, yAxis_max );
	c1->SetTickx(c1_orig->GetTickx());
	c1->SetTicky(c1_orig->GetTicky());
	c1->SetGridx(c1_orig->GetGridx());
	c1->SetGridy(c1_orig->GetGridy());
	c1->SetLogx(c1_orig->GetLogx());
	c1->SetLogy(c1_orig->GetLogy());
	c1->SetLogz(c1_orig->GetLogz());
	c1->SetBorderSize(c1_orig->GetBorderSize());
	c1->SetBorderMode(c1_orig->GetBorderMode());

	TString h_name_new = Form("%s", h_name.Data());
	TString h_title = Form("%s", title.Data());

	TH1D *h_base = new TH1D(h_name_new.Data(), h_title.Data(), f3_orig->GetNpx(), xMin, xMax);
	h_base->SetLineColor(0);
	h_base->SetStats(0);
	h_base->GetYaxis()->SetRangeUser(yAxis_min, yAxis_max);
	h_base->GetXaxis()->SetTitle(x_axis_title.Data());
	h_base->GetYaxis()->SetTitle(y_axis_title.Data());
	h_base->GetXaxis()->SetTitleOffset( h_orig_base->GetXaxis()->GetTitleOffset() );
	h_base->GetYaxis()->SetTitleOffset( h_orig_base->GetYaxis()->GetTitleOffset() );
	h_base->Draw();

	TF1 *f1 = (TF1*)f1_orig->Clone("f1");
	f1->SetTitle(title.Data());

	TF1 *f2 = (TF1*)f2_orig->Clone("f2");
	f2->SetTitle(title.Data());

	TF1 *f3 = (TF1*)f3_orig->Clone("f3");
	f3->SetTitle(title.Data());

	f1->SetFillColorAlpha(kRed - 10, 0.25);
	f1->SetFillStyle(1001);
	shade_btw_TF1s_w_3rd(c1, f1,f2,f3, yAxis_min, yAxis_max, f1_orig, f2_orig, f3_orig);
}

#endif