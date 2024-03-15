#ifndef PLOT_PIMPING_H
#define PLOT_PIMPING_H

struct tge_dataPoint{
	double x;
	double y;
	double y_error;
	double y_err_high;
	double y_err_low;
};

void plotShadedRegion(TGraphErrors *graph) {
    // Get the number of points in the graph
    const Int_t nPoints = graph->GetN();

    // Create a TGraphAsymmErrors with asymmetric errors
    TGraphAsymmErrors *asymmGraph = new TGraphAsymmErrors(nPoints);
    asymmGraph->SetName(graph->GetName());
    for (Int_t i = 0; i < nPoints; ++i) {
        Double_t x, y, ex, ey;
        graph->GetPoint(i, x, y);
        ex = graph->GetErrorX(i);
        ey = graph->GetErrorY(i);

        // Set asymmetric errors
        asymmGraph->SetPoint(i, x, y);
        asymmGraph->SetPointError(i, ex, ex, ey, ey);
    }

    // Create a TPolyLine for the shaded region
    TPolyLine *poly = new TPolyLine(2 * nPoints - 1);  // Use -1 to avoid connecting last to first point

    // Set points for the TPolyLine
    for (Int_t i = 0; i < nPoints - 1; ++i) {
        Double_t x, y, exLow, exHigh, eyLow, eyHigh;
        asymmGraph->GetPoint(i, x, y);
        exLow = asymmGraph->GetErrorXlow(i);
        exHigh = asymmGraph->GetErrorXhigh(i);
        eyLow = asymmGraph->GetErrorYlow(i);
        eyHigh = asymmGraph->GetErrorYhigh(i);

        // Set points for the TPolyLine
        poly->SetPoint(i, x + exLow, y - eyLow);
    }

    // Set points for the second segment of TPolyLine
    for (Int_t i = nPoints - 2; i >= 0; --i) {
        Double_t x, y, exHigh, eyHigh;
        asymmGraph->GetPoint(i, x, y);
        exHigh = asymmGraph->GetErrorXhigh(i);
        eyHigh = asymmGraph->GetErrorYhigh(i);

        poly->SetPoint(nPoints - 1 + (nPoints - 2 - i), x + exHigh, y + eyHigh);
    }

    // Create a canvas and draw the graph and shaded region
    TCanvas *canvas = new TCanvas("canvas", "Shaded Region Example", 800, 600);
    asymmGraph->Draw("AP");
    poly->SetFillColorAlpha(kBlue, 0.3);
    poly->Draw("F");

    // Add more customization or additional elements as needed
    canvas->Update();
}

TGraph *y_err_TGraphErrors(TGraphErrors &tge, TString err_select = ""){
	vector<tge_dataPoint> dataVector;
	Int_t nPoints = tge.GetN();
	vector<double> x_vec, y_vec, y_error_vec, y_err_low_vec, y_err_high_vec;
	x_vec.clear();
	y_vec.clear();
	y_error_vec.clear();
	y_err_low_vec.clear();
	y_err_high_vec.clear();

	dataVector.resize(nPoints);

	for( Int_t i = 0; i < nPoints; i++ ){
		double x = tge.GetX()[i];
		double y = tge.GetY()[i];
		double error = tge.GetEY()[i];
		double y_err_high = y + error;
		double y_err_low = y - error;

		x_vec.push_back(x);
		y_vec.push_back(y);
		y_error_vec.push_back(error);
		y_err_low_vec.push_back(y_err_low);
		y_err_high_vec.push_back(y_err_high);
	}

	TGraph *tge_err_highLow, *combinedGraph;
	if( err_select == "low" ){
		tge_err_highLow = new TGraph(nPoints, &x_vec[0], &y_err_low_vec[0]);
	} 
	else if( err_select == "high" ){
		tge_err_highLow = new TGraph(nPoints, &x_vec[0], &y_err_high_vec[0]);
	}
	else if( err_select == "combined" ){
		TGraph *tg1 = new TGraph(nPoints, &x_vec[0], &y_err_low_vec[0]);
		TGraph *tg2 = new TGraph(nPoints, &x_vec[0], &y_err_high_vec[0]);

		tge_err_highLow = new TGraph( tg1->GetN() + tg2->GetN());

		for( Int_t i = 0; i < tg1->GetN(); i++){
			tge_err_highLow->SetPoint(i, tg1->GetX()[i], tg1->GetY()[i]);
		}
		for( Int_t i = 0; i < tg2->GetN(); i++){
			tge_err_highLow->SetPoint(i, tg2->GetX()[i], tg2->GetY()[i]);
		}
		tge_err_highLow->SetPoint(tg1->GetN() + tg2->GetN(), tge_err_highLow->GetX()[0], tge_err_highLow->GetY()[0]);	
	}

	else{
		tge_err_highLow = new TGraph(nPoints, &x_vec[0], &y_vec[0]);
	}

	return tge_err_highLow;

}

void shade_between_TGraphs( TGraph *tg1, TGraph *tg2 ){
	TGraphPainter *painter = new TGraphPainter();

	painter->SetDrawOption("A3");
	// painter->SetFillColor(kBlue - 10);

	TGraph *combinedGraph = new TGraph( tg1->GetN() + tg2->GetN());

	for( Int_t i = 0; i < tg1->GetN()-1; i++){
		combinedGraph->SetPoint(i, tg1->GetX()[i], tg1->GetY()[i]);
	}
	for( Int_t i = 0; i < tg2->GetN()-1; i++){
		combinedGraph->SetPoint(tg1->GetN() + i, tg2->GetX()[i], tg2->GetY()[i]);
	}
	combinedGraph->Draw("A3");
	tg1->Draw("P+same");
	tg2->Draw("P+same");
}
void shade_btw_TGs_w_3rd(TCanvas *c1, TGraph *tg1, TGraph *tg2, TGraph *tg3, double y_Axis_min, double y_Axis_max) {
	//shade the area between f1 and f2 and draw f3 on top

	//create a TGraph to store the function values
	//shaded area is the fill/color/style of f1
	TGraph *gr = new TGraph();
	gr->SetFillColor(tg1->GetFillColor());
	gr->SetFillStyle(tg1->GetFillStyle());
	tg3->Draw("lsame");
	c1->Update();
	//get picture range
	Int_t npx = tg3->GetN();
	Double_t xmin = tg3->GetX()[0];
	Double_t xmax = tg3->GetX()[npx-1];
	Double_t ymin = y_Axis_min;
	Double_t ymax = y_Axis_max;

	//process first function

	Int_t npoints=0;
	Double_t dx = (xmax-xmin)/npx;
	for( int x = 0; x < npx; x++){
		Double_t y = tg1->GetY()[x];
		if (y < ymin) y = ymin;
		if (y > ymax) y = ymax;
		gr->SetPoint(npoints,x,y);
		npoints++;	
	}
	//process second function
	for( int x = 0; x < npx; x++){
		Double_t y = tg2->GetY()[x];
		if (y < ymin) y = ymin;
		if (y > ymax) y = ymax;
		gr->SetPoint(npoints,x,y);
		npoints++;
	}
	gr->Draw("fsame"); //draw graph with fill area option
	tg3->Draw("lsame"); //superimpose function
}
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

void fshade_btw_TF1s_w_3rd_SAME(TCanvas *c1, TH1D *h_orig_base, TF1 *f1_orig, TF1 *f2_orig, TF1 *f3_orig, double yAxis_min, double yAxis_max, int line_color = 0, int shade_color = 632) {
	TString h_name = h_orig_base->GetName();
	TString title = h_orig_base->GetTitle();
	TString x_axis_title = h_orig_base->GetXaxis()->GetTitle();
	TString y_axis_title = h_orig_base->GetYaxis()->GetTitle();

	TString canvas_name = Form("%s", h_name.Data());
	canvas_name.ReplaceAll("h_", "c_");

	TString canvas_title = Form("%s", title.Data());

	double xMin = f3_orig->GetXmin();
	double xMax = f3_orig->GetXmax();

	c1->cd();
	c1->Range(xMin, yAxis_min, xMax, yAxis_max );
	c1->SetTickx(c1->GetTickx());
	c1->SetTicky(c1->GetTicky());
	c1->SetGridx(c1->GetGridx());
	c1->SetGridy(c1->GetGridy());
	c1->SetLogx(c1->GetLogx());
	c1->SetLogy(c1->GetLogy());
	c1->SetLogz(c1->GetLogz());
	c1->SetBorderSize(c1->GetBorderSize());
	c1->SetBorderMode(c1->GetBorderMode());

	TString h_name_new = Form("%s", h_name.Data());
	TString h_title = Form("%s", title.Data());

	TH1D *h_base = new TH1D(h_name_new.Data(), h_title.Data(), f3_orig->GetNpx(), xMin, xMax);
	h_base->SetLineColor(line_color);
	h_base->SetStats(0);
	h_base->GetYaxis()->SetRangeUser(yAxis_min, yAxis_max);
	h_base->GetXaxis()->SetTitle(x_axis_title.Data());
	h_base->GetYaxis()->SetTitle(y_axis_title.Data());
	h_base->GetXaxis()->SetTitleOffset( h_orig_base->GetXaxis()->GetTitleOffset() );
	h_base->GetYaxis()->SetTitleOffset( h_orig_base->GetYaxis()->GetTitleOffset() );
	h_base->Draw("same");

	TF1 *f1 = (TF1*)f1_orig->Clone("f1");
	f1->SetTitle(title.Data());

	TF1 *f2 = (TF1*)f2_orig->Clone("f2");
	f2->SetTitle(title.Data());

	TF1 *f3 = (TF1*)f3_orig->Clone("f3");
	f3->SetTitle(title.Data());

	f1->SetFillColorAlpha(shade_color, 0.25);
	f1->SetFillStyle(1001);
	shade_btw_TF1s_w_3rd(c1, f1,f2,f3, yAxis_min, yAxis_max, f1_orig, f2_orig, f3_orig);
}

void drawSquareOnHistogram(TH2D* hist, Double_t x1, Double_t y1, Double_t x2, Double_t y2, int color = 2, TCanvas* canvas = nullptr) {
    // If canvas is not provided, create a new one
    if (!canvas)
        canvas = new TCanvas("canvas", "2D Histogram with Square", 800, 600);

    // Draw the histogram on the canvas
    hist->Draw("COLZ");

    // Define the square's edges
    Double_t x[5] = {x1, x2, x2, x1, x1};
    Double_t y[5] = {y1, y1, y2, y2, y1};

    // Create a polyline to draw the square
    TPolyLine* square = new TPolyLine(5, x, y);
    square->SetLineColor(color);
    square->SetLineWidth(2);
    square->Draw("same");

    // Update the canvas
    canvas->Update();
}

void drawSquareOnCanvas(TCanvas* canvas, Double_t x1, Double_t y1, Double_t x2, Double_t y2, int color = 2) {
    // If canvas is not provided, create a new one

    // Define the square's edges
    Double_t x[5] = {x1, x2, x2, x1, x1};
    Double_t y[5] = {y1, y1, y2, y2, y1};

    // Create a polyline to draw the square
    TPolyLine* square = new TPolyLine(5, x, y);
    square->SetLineColor(color);
    square->SetLineWidth(2);
    square->Draw("same");

    // Update the canvas
    canvas->Update();
}


#endif