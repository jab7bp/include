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

void build_and_set_branch_int(TChain *tChain, TString BranchVariable, Int_t &StoreVariable ){

	tChain->SetBranchStatus( BranchVariable.Data(), 1 );
	
	tChain->SetBranchAddress( BranchVariable.Data(), &StoreVariable );		
	// if( !Reference ){
	// 	tChain->SetBranchAddress( BranchName.Data(), StoreVariable );		
	// }

	cout << BranchVariable.Data() << " ";
}

void build_and_set_branch_double(TChain *tChain, TString BranchVariable, Double_t &StoreVariable ){

	tChain->SetBranchStatus( BranchVariable.Data(), 1 );
	
	tChain->SetBranchAddress( BranchVariable.Data(), &StoreVariable );		
	// if( !Reference ){
	// 	tChain->SetBranchAddress( BranchName.Data(), StoreVariable );		
	// }

	cout << BranchVariable.Data() << " ";
}

void build_and_set_branch_array(TChain *tChain, TString BranchVariable, Double_t *StoreVariable ){

	tChain->SetBranchStatus( BranchVariable.Data(), 1 );
	
	tChain->SetBranchAddress( BranchVariable.Data(), StoreVariable );		
	// if( !Reference ){
	// 	tChain->SetBranchAddress( BranchName.Data(), StoreVariable );		
	// }

	cout << BranchVariable.Data() << " ";
}


void FillVectorByFilesWithPattern(const char* folderPath, const TString& pattern, std::vector<TString>& matchingFileNames) {
    TSystemDirectory dir(folderPath, folderPath);

    TList* files = dir.GetListOfFiles();
    if (!files) {
        std::cerr << "Error: Unable to retrieve file list" << std::endl;
        return; // Return without modifying the vector
    }

    TIter next(files);
    TSystemFile* file;

    while ((file = (TSystemFile*)next())) {
        TString fileName(file->GetName());
        if (fileName.Contains(pattern)) {
        	// cout << "Found file: " << fileName.Data() << endl;
            matchingFileNames.push_back(fileName.Data());
        }
    }
}

void FillVectorByFilesAndDirWithPattern(const char* folderPath, const TString& pattern, std::vector<TString>& matchingFileNames) {
    TSystemDirectory dir(folderPath, folderPath);
    TString fileDir = Form("%s", folderPath);

    TList* files = dir.GetListOfFiles();
    if (!files) {
        std::cerr << "Error: Unable to retrieve file list" << std::endl;
        return; // Return without modifying the vector
    }

    TIter next(files);
    TSystemFile* file;

    while ((file = (TSystemFile*)next())) {
        TString fileName(file->GetName());
        if (fileName.Contains(pattern)) {
            // cout << "Found file: " << fileName.Data() << endl;
            matchingFileNames.push_back(Form("%s/%s", fileDir.Data(), fileName.Data() ));
        }
    }
}

bool CompareStringsWithNumbers(TString s1, TString s2)
{
 
    // If size of numeric strings
    // are same the put lowest value
    // first
    if (s1.Length() == s2.Length()) {
        return s1 < s2;
    }
 
    // If size is not same put the
    // numeric string with less
    // number of digits first
    else {
        return s1.Length() < s2.Length();
    }
}

void NormalizeHistogramToUnityIntegralSetBins(TH1D* histogram) {
    // Get the integral of the histogram
    Double_t integral = histogram->Integral();

    // Check if the integral is non-zero to avoid division by zero
    if (integral > 0.0) {
        // Normalize each bin content by the integral
        for (Int_t bin = 1; bin <= histogram->GetNbinsX(); ++bin) {
            histogram->SetBinContent(bin, histogram->GetBinContent(bin) / integral);
        }
    } else {
        std::cerr << "Warning: Histogram integral is zero. Cannot normalize." << std::endl;
    }
}

void NormalizeHistogramToUnityIntegralScaling(TH1D* histogram) {
    // Get the integral of the histogram
    Double_t integral = histogram->Integral();

    // Check if the integral is non-zero to avoid division by zero
    if (integral > 0.0) {
        // Normalize each bin content by the integral
        histogram->Scale(1.0/integral);
    } else {
        std::cerr << "Warning: Histogram integral is zero. Cannot normalize." << std::endl;
    }
}

double searchSimcHistFile( const TString &searchPattern, const TString &filename ){
    std::ifstream inputFile(filename.Data()); // Open the file for reading

    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        std::cerr << filename.Data() << std::endl;
        return -99;
    }

    TString line;
    double foundValue = -1.0; // Default value if not found, using a double
    TString previousLine; // Store the previous line for exceptions
    while (line.ReadLine(inputFile)) {
        // Search for lines containing the specified search pattern
        if (line.Contains(searchPattern)) {
            // Exception for "luminosity" and other patterns with scientific notation
            if (line.EndsWith("ub^-1") || line.EndsWith("geV^2")) {
                // Extract the numerical value using TString::Tokenize
                TString delim("=");
                TObjArray* tokens = line.Tokenize(delim);
                if (tokens->GetEntries() >= 2) {
                    TObjString* valueStr = dynamic_cast<TObjString*>(tokens->At(1));
                    if (valueStr) {
                        TString valueString = valueStr->String();
                        // Convert the value to a double
                        foundValue = valueString.Atof();
                        delete tokens; // Clean up the tokens
                        return foundValue;
                    }
                }
                delete tokens; // Clean up the tokens if no valid value was found
            } else {
                // For other patterns, try to parse as a double
                TString valueString;
                TString delim("=");
                TObjArray* tokens = line.Tokenize(delim);
                if (tokens->GetEntries() >= 2) {
                    TObjString* valueStr = dynamic_cast<TObjString*>(tokens->At(1));
                    if (valueStr) {
                        valueString = valueStr->String();
                    }
                }
                // Attempt to convert to a double
                if (valueString.IsFloat()) {
                    foundValue = valueString.Atof();
                    delete tokens; // Clean up the tokens
                    return foundValue;
                }
                delete tokens; // Clean up the tokens if no valid value was found
            }
        }
        // Store the previous line for exceptions
        previousLine = line;
    }

    // Handle the "luminosity" and other exceptions when no appropriate line is found
    if ((searchPattern == "luminosity" && previousLine.EndsWith("ub^-1")) ||
        (searchPattern == "genvol" )) {
        // Extract the numerical value using TString::Tokenize
        TString delim("=");
        TObjArray* tokens = previousLine.Tokenize(delim);
        if (tokens->GetEntries() >= 2) {
            TObjString* valueStr = dynamic_cast<TObjString*>(tokens->At(1));
            if (valueStr) {
                TString valueString = valueStr->String();
                // Convert the value to a double
                foundValue = valueString.Atof();
            }
        }
        delete tokens; // Clean up the tokens if no valid value was found
    }

    inputFile.close(); // Close the file

    return foundValue; // Return the default value if not found
}


//Function to find the "Best Cluster" for hcal_e given various criteria

//For these we need to feed it the vectors to pick from:
//Clusters which pass the coin timing
//Energy and timing value/diff
//Corresponding indexes
void bestClusterSelectionHcalE_CoinPassAndMaxE(
    const std::vector<Double_t>& energy_vec,
    const std::vector<Int_t>& index_vec,
    std::vector<Double_t>& hcal_cluster_energy_and_index
) {
    if (energy_vec.empty() || index_vec.empty()) {
        std::cerr << "Error: Input vectors are empty." << std::endl;
        return;
    }

    // Find the index of the maximum value in energy_vec
    double maxEnergy = energy_vec[0];
    int maxIndex = 0;

    for (size_t i = 1; i < energy_vec.size(); ++i) {
        if (energy_vec[i] > maxEnergy) {
            maxEnergy = energy_vec[i];
            maxIndex = i;
        }
    }

    // Retrieve the value of index_vec at the maxIndex and store it as a double
    double maxIndexValue = static_cast<double>(index_vec[maxIndex]);

    // Populate the result vector
    hcal_cluster_energy_and_index.clear();
    hcal_cluster_energy_and_index.push_back(maxIndexValue);
    hcal_cluster_energy_and_index.push_back(maxEnergy);
    hcal_cluster_energy_and_index.push_back(maxIndexValue);
}

void bestClusterSelectionHcalE_CoinTimingMaxE_Scoring(
    const std::vector<Double_t>& energy_vec,
    const std::vector<Int_t>& index_vec,
    const std::vector<Double_t>& timing_vec,
    std::vector<Double_t>& hcal_cluster_energy_and_index
) {
    if (energy_vec.empty() || index_vec.empty() || timing_vec.empty()) {
        std::cerr << "Error: Input vectors are empty." << std::endl;
        sleep(10);
        exit(1);
        return;
    }

    // Initialize variables to store the best combination
    double bestScore = -1.0;  // Initialize to a negative value
    double bestEnergy = 0.0;
    double bestTiming = 0.0;
    int bestIndex = 0;
    vector<double> scores = {};
    int originalIndex = -1;

    // Evaluate each combination of energy and timing
    for (size_t i = 0; i < energy_vec.size(); ++i) {
        double score = energy_vec[i] - timing_vec[i];
        scores.push_back( score );

        // If the current combination has a better score, update the best values
        if (score > bestScore || (score == bestScore && timing_vec[i] < bestTiming)) {
            bestScore = score;
            bestEnergy = energy_vec[i];
            bestTiming = timing_vec[i];
            bestIndex = index_vec[i];
            originalIndex = i;
        }
    }

    // Populate the result vector
    hcal_cluster_energy_and_index.clear();
    hcal_cluster_energy_and_index.push_back(originalIndex);
    hcal_cluster_energy_and_index.push_back(bestEnergy);
    hcal_cluster_energy_and_index.push_back(static_cast<double>(bestIndex));
    hcal_cluster_energy_and_index.push_back(bestTiming);
    hcal_cluster_energy_and_index.push_back(bestScore);
    // for( size_t i = 0; i < scores.size(); i++ ){
    //     cout << scores[i] << " ";
    // }
}

void normalizeVector(std::vector<Double_t>& values) {
    if (values.empty()) {
        return; // Avoid division by zero if the vector is empty
    }

    // Find the minimum and maximum values in the vector
    Double_t min_val = *std::min_element(values.begin(), values.end());
    Double_t max_val = *std::max_element(values.begin(), values.end());

    // Normalize the values to the range [0, 1]
    for (Double_t& value : values) {
        value = (value - min_val) / (max_val - min_val);
    }
}

void bestClusterSelectionHcalE_CoinTimingMaxE_NormalizedScoring(
    std::vector<Double_t>& energy_vec,
    const std::vector<Int_t>& index_vec,
    std::vector<Double_t>& timing_vec,
    std::vector<Double_t>& hcal_cluster_energy_and_index
) {
    if (energy_vec.empty() || index_vec.empty() || timing_vec.empty()) {
        std::cerr << "Error: Input vectors are empty." << std::endl;
        sleep(10);
        exit(1);
        return;
    }

    // Initialize variables to store the best combination
    double bestScore = -1.0;  // Initialize to a negative value
    double bestEnergy = 0.0;
    double bestTiming = 0.0;
    int bestIndex = 0;
    vector<double> scores = {};
    int originalIndex = -1;

    // Normalize the energy and timing values
    normalizeVector(energy_vec);
    normalizeVector(timing_vec);

    // Evaluate each combination of energy and timing
    for (size_t i = 0; i < energy_vec.size(); ++i) {
        Double_t score = energy_vec[i] + timing_vec[i];
        scores.push_back( score );

       // If the current combination has a better score, update the best values
        if (score > bestScore) {
            bestScore = score;
            bestIndex = i;
            bestEnergy = energy_vec[i];
            bestTiming = timing_vec[i];
            originalIndex = i;
        }
    }

    // Populate the result vector
    hcal_cluster_energy_and_index.clear();
    hcal_cluster_energy_and_index.push_back(originalIndex);
    hcal_cluster_energy_and_index.push_back(bestEnergy);
    hcal_cluster_energy_and_index.push_back(static_cast<double>(bestIndex));
    hcal_cluster_energy_and_index.push_back(bestTiming);
    hcal_cluster_energy_and_index.push_back(bestScore);
    // for( size_t i = 0; i < scores.size(); i++ ){
    //     cout << scores[i] << " ";
    // }
}

Double_t euclideanDistance(Double_t x, Double_t y) {
    return std::sqrt(x * x + y * y);
}

void bestClusterSelectionHcalE_CoinTimingMaxE_EuclideanScoring(
    const std::vector<Double_t>& energy_vec,
    const std::vector<Int_t>& index_vec,
    const std::vector<Double_t>& timing_vec,
    std::vector<Double_t>& hcal_cluster_energy_and_index
) {
    if (energy_vec.empty() || index_vec.empty() || timing_vec.empty()) {
        std::cerr << "Error: Input vectors are empty." << std::endl;
        return;
    }

    // Find the best combination based on Euclidean Distance
    int bestIndex = -1;
    Double_t bestDistance = std::numeric_limits<Double_t>::max(); // Initialize to a large value

    // Initialize variables to store the best combination
    double bestEnergy = 0.0;
    double bestTiming = 0.0;
    int originalIndex = -1;

    for (size_t i = 0; i < energy_vec.size(); ++i) {
        Double_t distance = euclideanDistance(energy_vec[i], timing_vec[i]);

        // If the current combination has a smaller distance, update the best values
        if (distance < bestDistance) {
            bestDistance = distance;
            bestIndex = i;
            bestEnergy = energy_vec[i];
            bestTiming = timing_vec[i];
            originalIndex = i;
        }
    }

    // Populate the result vector with bestIndex first
    hcal_cluster_energy_and_index.clear();
    hcal_cluster_energy_and_index.push_back(originalIndex);
    hcal_cluster_energy_and_index.push_back(bestEnergy);
    hcal_cluster_energy_and_index.push_back(static_cast<double>(bestIndex));
    hcal_cluster_energy_and_index.push_back(bestTiming);

}

void shiftHistogram(TH1D* hist, Double_t x_shift) {
    // Check if the histogram is valid
    if (!hist) {
        std::cerr << "Error: Invalid histogram." << std::endl;
        return;
    }

    // Create a new histogram with the same properties as the original
    TH1D* shiftedHist = dynamic_cast<TH1D*>(hist->Clone());

    // Loop over all bins and shift the x-axis values
    for (Int_t bin = 1; bin <= hist->GetNbinsX(); ++bin) {
        Double_t binCenter = hist->GetBinCenter(bin);
        Double_t shiftedBinCenter = binCenter + x_shift;
        Double_t binContent = hist->GetBinContent(bin);

        // Find the corresponding bin for the shifted bin center
        Int_t shiftedBin = shiftedHist->FindBin(shiftedBinCenter);

        // Set the content in the shifted histogram
        shiftedHist->SetBinContent(shiftedBin, binContent);
    }

    // Copy the contents of the shifted histogram back to the original histogram
    hist->Reset();  // Clear the original histogram
    hist->Add(shiftedHist);  // Add the shifted histogram contents to the original

    // Update the histogram
    hist->SetDirectory(0);  // Disassociate from current ROOT directory
    // hist->SetXTitle(Form("%s + %.2f", hist->GetXaxis()->GetTitle(), x_shift));

    // Clean up the dynamically allocated histogram
    delete shiftedHist;
}

void addHistogramToRootFile(const char* fileName, TH1* histogram, const char* directoryName) {
    // Open the ROOT file in update mode, creating it if it doesn't exist
    TFile* rootFile = new TFile(fileName, "UPDATE");

    // Check if the file is open and writable
    if (!rootFile || !rootFile->IsOpen() || rootFile->IsZombie()) {
        std::cerr << "Error: Unable to open or create the ROOT file." << std::endl;
        return;
    }

    // Check if the directory exists, and create it if necessary
    TDirectory* dir = rootFile->GetDirectory(directoryName);
    if (!dir) {
        dir = rootFile->mkdir(directoryName);
    }

    // Switch to the directory
    dir->cd();

    // Write the histogram to the file
    histogram->Write();

    // Close the file
    rootFile->Close();

    // Clean up
    delete rootFile;
}

double IntegralPostiveValsOnlyInRangeX(TH1* histogram, double x_min, double x_max) {
    // Get the histogram's x-axis
    TAxis* xAxis = histogram->GetXaxis();

    int binMin = xAxis->FindBin(x_min);
    int binMax = xAxis->FindBin(x_max);

    // Get the number of bins
    int numBins = xAxis->GetNbins();

    // Initialize the sum
    double sum = 0.0;

    // Loop over the bins and sum the positive values
    for (int bin = binMin; bin <= binMax; ++bin) {
        double binContent = histogram->GetBinContent(bin);
        if (binContent > 0) {
            sum += histogram->GetBinContent(bin);
        }
    }

    return sum;
}


#endif