#ifndef utility_functions
#define utility_functions

struct Array1DValueWithIndex {
	double ArrayValue;
	int index;
};

struct Array1DScoredWithIndex {
	int clusterIndex;
	int energyIndex;
	int timingIndex;
	int dxdyIndex;
};

void PrintScoredArray( Array1DScoredWithIndex &PrintStruct ){
	cout << "{" << PrintStruct.clusterIndex << ", " << PrintStruct.energyIndex << ", " << PrintStruct.timingIndex << ", " << PrintStruct.dxdyIndex << "}" << endl;
}

bool CompareArrayValueWithIndex( const Array1DValueWithIndex &a, const Array1DValueWithIndex &b ){
	return a.ArrayValue > b.ArrayValue; //Sort in descending order
}

//Fucntion to create a new, sorted array while also keeping track of the original indices before the sorting
void SortArrayWithIndices(double *originalArray, int arraySize, Array1DValueWithIndex *sortedArray ){

	for( int i = 0; i < arraySize; i++ ){
		sortedArray[i].ArrayValue = originalArray[i];
		sortedArray[i].index = i;
	}
	sort( sortedArray, sortedArray + arraySize, CompareArrayValueWithIndex );

}

//Function to sort an array while also keeping track of the original indices before the sorting. 
//MODIFIES ORIGINAL ARRAY
void SortAndModifyArrayWithIndices( vector<double> &values, vector<int> indices ){

	//Create a vector of pairs where each pair contains the original value;
	vector<pair<double, int>> valueIndexPairs;

	for( size_t i = 0; i < values.size(); i++ ){
		valueIndexPairs.push_back(make_pair(values[i], i));
	}

	//Sort the vector in pairs in descending order based on the values
	sort( valueIndexPairs.rbegin(), valueIndexPairs.rend() );

	//Update the original arrays with the sorted values and corresponding indices
	for( size_t i = 0; i < values.size(); i++ ){
		values[i] = valueIndexPairs[i].first;
		indices[i] = valueIndexPairs[i].second;
	}

}

int FindIndexOfMaxElement( const double array[], int size ){
	if (size <= 0) {
        // Handle empty array or invalid size
        return -1; // Return a sentinel value or an error code
    }

    int maxIndex = 0; // Initialize maxIndex to the index of the first element

    for (int i = 1; i < size; ++i) {
        if (array[i] > array[maxIndex]) {
            maxIndex = i; // Update maxIndex if a larger element is found
        }
    }

    return maxIndex;

}

TH1D *make_residual_histo(TString identifier, TH1D *histo_1, TH1D *histo_2, bool match_bins = true, bool match_x_range = false ){

	TString histo_1_name = histo_1->GetName();
	Int_t Nbins_histo_1 = histo_1->GetNbinsX();
	Double_t MinX_histo_1 = histo_1->GetXaxis()->GetXmin();
	Double_t MaxX_histo_1 = histo_1->GetXaxis()->GetXmax();

	TString histo_2_name = histo_2->GetName();
	Int_t Nbins_histo_2 = histo_2->GetNbinsX();
	Double_t MinX_histo_2 = histo_2->GetXaxis()->GetXmin();
	Double_t MaxX_histo_2 = histo_2->GetXaxis()->GetXmax();

	TH1D *resid_histo = new TH1D( "resid_histo", Form("Residual Histogram - %s: %s and %s", identifier.Data(), histo_1_name.Data(), histo_2_name.Data()), Nbins_histo_1, MinX_histo_1, MaxX_histo_1);

	if( match_bins && (Nbins_histo_1 != Nbins_histo_2) ){
		cout << "------------------------------------------------------" << endl;
		cout << "------------------------------------------------------" << endl;
		cout << "---- Histograms need to have matching bin numbers ----" << endl;
		cout << "     Histo_1: " << Nbins_histo_1 << ", Histo_2: " << Nbins_histo_2 << endl;
		cout << "------------------------------------------------------" << endl;
		cout << "------------------------------------------------------" << endl;
		return resid_histo;
	}

	if( match_x_range && (MinX_histo_1 != MinX_histo_2) && (MaxX_histo_1 != MaxX_histo_2 )){
		cout << "------------------------------------------------------" << endl;
		cout << "------------------------------------------------------" << endl;
		cout << "----- Histograms need to have matching x ranges  -----" << endl;
		cout << "Min X - Histo_1: " << MinX_histo_1 << ", Histo_2: " << MinX_histo_2 << endl;
		cout << "Max X - Histo_1: " << MaxX_histo_1 << ", Histo_2: " << MaxX_histo_2 << endl;
		cout << "------------------------------------------------------" << endl;
		cout << "------------------------------------------------------" << endl;
		return resid_histo;		
	}

	for( int bin = 0; bin < Nbins_histo_1; bin++ ){
		resid_histo->SetBinContent( bin, histo_1->GetBinContent(bin) - histo_2->GetBinContent(bin) );
	}

	return resid_histo;
}

#endif