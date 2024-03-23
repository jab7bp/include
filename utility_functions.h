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


bool CompareFilenames(const TString& filename1, const TString& filename2) {
    int val1_start = filename1.Last('_')+1;
    int val1_end = filename1.Last('.') - 1;
    int val1_length = val1_end - val1_start;

    int val2_start = filename2.Last('_')+1;
    int val2_end = filename2.Last('.') - 1;
    int val2_length = val2_end - val2_start;
 
    TString digits1 = filename1(val1_start, val1_length);
    TString digits2 = filename2(val2_start, val2_length);

    // TString digits1 = filename1(filename1.Length() - 8, 4);
    // TString digits2 = filename2(filename2.Length() - 8, 4);

    return digits1.Atoi() < digits2.Atoi();
}

// bool compareByNumber(const TString& str1, const TString& str2) {
//     // Find the position of the last underscore and the dot before the extension
//     int val1_start = str1.Last('_');
//     int val1_end = str1.Last('.') - 1;
//     int val1_length = val1_end - val1_start;

//     int val2_start = str2.Last('_');
//     int val2_end = str2.Last('.') - 1;
//     int val2_length = val2_end - val2_start;
 
//     TString digits1 = str1(val1_start, val1_length);
//     TString digits2 = str2(val2_start, val2_length);

//     int num1 = digits1.Atoi();
//     int num2 = digits2.Atoi();
    
//     // Compare the numbers
//     return num1 < num2;
// }

// bool compareByStringNumber(const TString& str1, const TString& str2) {
//     // Find the position of the last underscore and the dot before the extension
//     Ssiz_t pos1 = str1.Last('_');
//     Ssiz_t pos2 = str2.Last('_');
//     Ssiz_t dotPos1 = str1.Last('.');
//     Ssiz_t dotPos2 = str2.Last('.');
    
//     // Extract the numbers from the substrings
//     TString numStr1 = str1(pos1 + 1, dotPos1 - pos1 - 1);
//     TString numStr2 = str2(pos2 + 1, dotPos2 - pos2 - 1);
    
//     // Compare the numbers as strings
//     return numStr1.CompareTo(numStr2) < 0;
// }

// Function to extract number from string
int extractNumber(const TString& str) {
    Ssiz_t pos = str.Last('_');
    Ssiz_t dotPos = str.Last('.');
    TString numStr = str(pos + 1, dotPos - pos - 1);
    return numStr.Atoi();
}

// Custom comparison function
bool compareStringsByLastNumber(const TString& str1, const TString& str2) {
    int num1 = extractNumber(str1);
    int num2 = extractNumber(str2);
    
    return num1 < num2;
}


void FillTwoVectorsByFilesInTwoDirsWithPatternAndSort(const TString& lookupPattern, const TString& dir1, const TString& dir2, std::vector<TString>& files1, std::vector<TString>& files2) {
    //Made for systematics analysis
    //Looks to pull a data and simc file that have the same systematics cut. 
    //The filenames are the same here but located in different locations.
    //The last four digits after the last underscore should match between each file
    //The last four digits correspond to some specific cut value or range.
    //Could be modified to accept two values separated by an underscore --> Represents a range or window of values

    //Will fill two pre-created vectors. 

    // Clear the vectors
    files1.clear();
    files2.clear();

    // Get the list of files in dir1
    TSystemDirectory sysDir1(dir1, dir1);
    // Get the list of files in dir2
    TSystemDirectory sysDir2(dir2, dir2);

    TList *filesList1 = sysDir1.GetListOfFiles();
    TList *filesList2 = sysDir2.GetListOfFiles();

    if (filesList1) {
        filesList1->Sort(&CompareFilenames); // Sort the list by last 4 digits
        TSystemFile *file;
        TString filename;
        TIter next(filesList1);
        while ((file = (TSystemFile*)next())) {
            filename = file->GetName();
            // cout << "Data filename: " << filename.Data() << endl;
            if (!file->IsDirectory() && filename.EndsWith(".root")) {
                TString lastFourDigits = filename(filename.Length() - 8, 4);
                TString partnerFilename = Form("%s/%s", dir2.Data(), filename.Data());
                // partnerFilename[partnerFilename.Length() - 8, 4] = lastFourDigits; // Replace last 4 digits
                if (gSystem->AccessPathName(partnerFilename) || !partnerFilename.Contains(lookupPattern))
                    continue;
                files1.push_back(dir1 + "/" + filename);
                files2.push_back(partnerFilename);
            }
        }
    }

    sort(files1.begin(), files1.end(), compareStringsByLastNumber);
    sort(files2.begin(), files2.end(), compareStringsByLastNumber);

    // Clean up
    delete filesList1;
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

bool CompareStringsByLastNumberAfterUnderscore(TString s1, TString s2){
    //First we get the positions of the last underscores for both strings
    Ssiz_t s1_lastUnderscorePos = s1.Last('_');
    Ssiz_t s2_lastUnderscorePos = s2.Last('_');

    TString s1_valueString = s1(s1_lastUnderscorePos + 1, s1.Last('.') - s1_lastUnderscorePos - 1 );
    TString s2_valueString = s2(s2_lastUnderscorePos + 1, s2.Last('.') - s2_lastUnderscorePos - 1 );
    
    //Convert those values to a double:
    double s1_extractedValue = s1_valueString.Atof();
    double s2_extractedValue = s2_valueString.Atof();

    //Compare....
    // If size of the value is the same (shouldn't be), put s1 first I guess...
    if( s1_extractedValue == s2_extractedValue ){
        return s1 < s2;
    }

    //otherwise, return the true lower value anyhow
    else{
        return s1_extractedValue < s2_extractedValue;
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

void shiftHistogram_by_kine( TH1D *histo_to_shift, int kine, bool systematics_mode = false ){

    if( systematics_mode ){
         if( kine == 4 ){
            shiftHistogram(histo_to_shift, 0.05);
        }

        if( kine == 8 ){
            shiftHistogram(histo_to_shift, .00);
        }

        if( kine == 9 ){
            shiftHistogram(histo_to_shift, -0.01);
        }   
    }
    else{

        if( kine == 4 ){
            shiftHistogram(histo_to_shift, 0.00);
        }

        if( kine == 8 ){
            shiftHistogram(histo_to_shift, 0.01);
        }

        if( kine == 9 ){
            shiftHistogram(histo_to_shift, -0.025);
        }   
    }


}


TH1D *make_shifted_histogram( TH1D *histo_to_shift, int kine ){

    TString histo_name = histo_to_shift->GetName();
    TString histo_title = histo_to_shift->GetTitle();

    TH1D *shifted_histogram = (TH1D*)histo_to_shift->Clone(histo_name.Data());
    shifted_histogram->SetTitle( histo_title.Data() );
    
    if( kine == 4 ){
        shiftHistogram(shifted_histogram, 0.04);
    }

    if( kine == 8 ){
        shiftHistogram(shifted_histogram, 0.03);
    } 

    return shifted_histogram;
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


TH1D* convertHistogramTo500Bins(const TH1D* sourceHistogram, const char* newHistogramName) {
    // Determine the number of bins in the source histogram
    int currentBins = sourceHistogram->GetNbinsX();
    
    // Determine the desired number of bins
    int targetBins = 500;

    // Create a new histogram with the desired number of bins
    TH1D* targetHistogram = new TH1D(newHistogramName, "", targetBins, sourceHistogram->GetXaxis()->GetBinLowEdge(1), sourceHistogram->GetXaxis()->GetBinUpEdge(currentBins));

    // Loop over the new bins and interpolate/average content from the original bins
    for (int i = 1; i <= targetBins; ++i) {
        double newBinCenter = targetHistogram->GetXaxis()->GetBinCenter(i);
        int bin = sourceHistogram->GetXaxis()->FindFixBin(newBinCenter);
        targetHistogram->SetBinContent(i, sourceHistogram->GetBinContent(bin));
    }

    return targetHistogram;
}


TH1D* createHistogramWithDoubleBinWidth(const TH1D* originalHist) {
    // Get information from the original histogram
    int originalBins = originalHist->GetNbinsX();
    double originalXmin = originalHist->GetXaxis()->GetXmin();
    double originalXmax = originalHist->GetXaxis()->GetXmax();

    // Create a new histogram with double the bin width
    int newBins = originalBins / 2;
    double newBinWidth = (originalXmax - originalXmin) / newBins;
    TH1D* newHist = new TH1D("newHist", "Histogram with Double Bin Width", newBins, originalXmin, originalXmax + newBinWidth);

    // Fill the new histogram with data from the original histogram
    for (int i = 1; i <= newBins; ++i) {
        double binContent = 0.0;
        double binError = 0.0;

        // Combine adjacent bins in the original histogram
        for (int j = 0; j < 2; ++j) {
            int originalBin = 2 * (i - 1) + 1 + j;
            binContent += originalHist->GetBinContent(originalBin);
            binError += originalHist->GetBinError(originalBin) * originalHist->GetBinError(originalBin);
        }

        binError = sqrt(binError);

        // Set the bin content and error in the new histogram
        newHist->SetBinContent(i, binContent);
        newHist->SetBinError(i, binError);
    }

    return newHist;
}

TH1D* createHistogramWithCustomBinWidth(const TH1D* originalHist, double binWidthMultiplier) {
    // Get information from the original histogram
    double original_maximum = originalHist->GetMaximum();

    int originalBins = originalHist->GetNbinsX();
    double originalXmin = originalHist->GetXaxis()->GetXmin();
    double originalXmax = originalHist->GetXaxis()->GetXmax();

    // Calculate the number of bins for the new histogram
    int newBins = originalBins / binWidthMultiplier;
    double newBinWidth = (originalXmax - originalXmin) / newBins;

    // Create a new histogram with custom bin width
    TH1D* newHist = new TH1D("newHist", "Histogram with Custom Bin Width", newBins, originalXmin, originalXmax + newBinWidth);

    // Fill the new histogram with data from the original histogram
    for (int i = 1; i <= newBins; ++i) {
        double binContent = 0.0;
        double binError = 0.0;

        // Combine adjacent bins in the original histogram
        for (int j = 0; j < binWidthMultiplier; ++j) {
            int originalBin = binWidthMultiplier * (i - 1) + 1 + j;
            binContent += originalHist->GetBinContent(originalBin);
            binError += originalHist->GetBinError(originalBin) * originalHist->GetBinError(originalBin);
        }

        binError = sqrt(0.03*newBins*binError);

        // Set the bin content and error in the new histogram
        newHist->SetBinContent(i, binContent);
        newHist->SetBinError(i, binError);
    }
    double new_maximum = newHist->GetMaximum();
    double normalize_scale_factor = original_maximum/new_maximum;

    newHist->Scale(normalize_scale_factor);

    return newHist;
}

int findTextFileLineNumberWithCustomString(const TString& filename, const TString& searchText) {
    std::ifstream inputFile(filename.Data());

    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open the file " << filename.Data() << std::endl;
        return -1; // Return -1 to indicate an error
    }

    TString line;
    int lineNumber = 0;

    while (inputFile.good()) {
        line.ReadLine(inputFile);
        lineNumber++;

        // Check if the line contains the search text
        if (line.Contains(searchText)) {
            inputFile.close();
            return lineNumber;
        }
    }

    inputFile.close();

    // Return -1 if the search text is not found
    return -1;
}

bool replaceTextInLine(const TString& filename, int lineNumber, const TString& newText) {
    std::ifstream inputFile(filename.Data());

    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open the file " << filename.Data() << std::endl;
        return false; // Return false to indicate an error
    }

    TString line;
    TString fileContent;
    int currentLineNumber = 0;

    while (inputFile.good()) {
        line.ReadLine(inputFile);
        currentLineNumber++;

        // Check if this is the line to replace
        if (currentLineNumber == lineNumber) {
            line = newText;
        }

        fileContent += line + "\n"; // Add a newline character after each line
    }

    inputFile.close();

    // Write the modified content back to the file
    std::ofstream outputFile(filename.Data());

    if (!outputFile.is_open()) {
        std::cerr << "Error: Unable to open the file " << filename.Data() << " for writing." << std::endl;
        return false; // Return false to indicate an error
    }

    outputFile << fileContent;
    outputFile.close();

    return true; // Return true if the replacement was successful
}

void processSystematicTextFile( const TString& filename, const TString& searchText, const TString& newText ){
    int line_number = findTextFileLineNumberWithCustomString(filename, searchText);
    if( line_number != -1 ){
        replaceTextInLine(filename, line_number, newText);
    }
    if( line_number == -1 ){
        std::ofstream outputFile(filename.Data(), std::ios::app );
        outputFile << newText.Data() << endl;
        outputFile.close();
    }
}

vector<double> ExtractNumericalToRangeFromFilename(TString filename) {
    TString search_anchor = "_to_";
    int length_of_vals = 4;

    double val1, val2;
    TString str1, str2;

    int anchor_index = filename.Index(search_anchor);

    str1 = filename(anchor_index - length_of_vals, length_of_vals);
    str2 = filename(anchor_index + search_anchor.Length(), length_of_vals);

    val1 = str1.Atof()/1000.0;
    val2 = str2.Atof()/1000.0;

    vector< double> values = {val1, val2};

    return values;

}

double ExtractWindowSizeFromFilename(TString filename){
    TString search_anchor = "windowSize";
    int length_of_vals = 4;

    int anchor_index = filename.Index(search_anchor);

    TString str = filename(anchor_index + search_anchor.Length(), length_of_vals);

    double val = str.Atof()/1000.0;
    return val;
}

double findMinWithErrors(const std::vector<double>& data, const std::vector<double>& errors) {
    double min_with_errors = data[0] - errors[0]; // Initialize with the first data point and its error bar
    for (size_t i = 1; i < data.size(); ++i) {
        double value_with_error = data[i] - errors[i];
        if (value_with_error < min_with_errors) {
            min_with_errors = value_with_error;
        }
    }
    return min_with_errors;
}

// Function to find maximum value considering the error bars
double findMaxWithErrors(const std::vector<double>& data, const std::vector<double>& errors) {
    double max_with_errors = data[0] + errors[0]; // Initialize with the first data point and its error bar
    for (size_t i = 1; i < data.size(); ++i) {
        double value_with_error = data[i] + errors[i];
        if (value_with_error > max_with_errors) {
            max_with_errors = value_with_error;
        }
    }
    return max_with_errors;
}

void sortTwoCorrelatedVectors(std::vector<double>& vec1, std::vector<double>& vec2) {
    // Create a vector of pairs where the first element is the value from vec1
    // and the second element is the index
    std::vector<std::pair<double, int>> pairs;
    for (int i = 0; i < vec1.size(); ++i) {
        pairs.push_back({vec1[i], i});
    }

    // Sort the vector of pairs based on the value (first element)
    std::sort(pairs.begin(), pairs.end());

    // Apply the same permutation to vec1 and vec2
    std::vector<double> sortedVec1;
    std::vector<int> sortedVec2;
    for (const auto& pair : pairs) {
        sortedVec1.push_back(pair.first);
        sortedVec2.push_back(vec2[pair.second]);
    }

    // Update the original vectors
    // vec1 = sortedVec1;
    // vec2 = sortedVec2;
    for(int i = 0; i < sortedVec1.size(); i++ ){
        vec1[i] = sortedVec1[i];
        vec2[i] = sortedVec2[i];
    }
}

double linear_interpolation( double x1, double y1, double x2, double y2, double x ){
    double interpolation = y1 + (x - x1)*( (y2 - y1)/(x2 - x11) );

    return interpolation;
}

std::vector<std::pair<double, double>>  GetHistogramXvalSlicesByNevents(TH2D* hist, double targetEvents, double x_min = NAN, double x_max = NAN) {
    // Let's allow for some custom start/finish x values:
    bool custom_x_bounds = false;

    int numBinsX = hist->GetNbinsX();
    int numBinsY = hist->GetNbinsY();

    int startBin = 1;
    int endBin = numBinsX;

    // if x_min and x_max are still NaN then do nothing, else, grab their corresponding bin values and use thouse as bounds
    if( (x_min != x_min) || (x_max != x_max) ){
        custom_x_bounds = false;
    }
    else{
        custom_x_bounds = true;
        startBin = hist->GetXaxis()->FindBin(x_min);
        endBin = hist->GetXaxis()->FindBin(x_max);
    }

    // Iterate over x-axis bins
    double totalEvents = 0;
    double sliceStart = hist->GetXaxis()->GetBinLowEdge(startBin);
    std::vector<std::pair<double, double>> slices; // Stores (start, end) pairs of x-slices

    for (int i = startBin; i <= endBin; ++i) {
        // Sum events in this x-bin
        for (int j = 1; j <= numBinsY; ++j) {
            totalEvents += hist->GetBinContent(i, j);
        }

        // Check if target events reached
        if (totalEvents >= targetEvents || i == numBinsX) {
            double sliceEnd = hist->GetXaxis()->GetBinUpEdge(i);

            // cout << "Found range: " << sliceStart << " -> " << sliceEnd << ", " << totalEvents << endl;
            // cout << "-------------------" << endl;
            // Store the slice
            slices.push_back(std::make_pair(sliceStart, sliceEnd));

            // Reset for next slice
            totalEvents = 0;
            sliceStart = sliceEnd;
        }
    }
    return slices;

}

std::vector<std::pair<double, double>>  GetHistogramBinSlicesByNevents(TH2D* hist, double targetEvents, double x_min = NAN, double x_max = NAN) {
    bool print = false;

    // Let's allow for some custom start/finish x values:
    bool custom_x_bounds = false;
    
    int numBinsX = hist->GetNbinsX();
    int numBinsY = hist->GetNbinsY();  

    int startBin = 1;
    int endBin = numBinsX;

    // if x_min and x_max are still NaN then do nothing, else, grab their corresponding bin values and use thouse as bounds
    if( (x_min != x_min) || (x_max != x_max) ){
        custom_x_bounds = false;
    }
    else{
        custom_x_bounds = true;
        startBin = hist->GetXaxis()->FindBin(x_min);
        endBin = hist->GetXaxis()->FindBin(x_max);
    }


    // Iterate over x-axis bins
    double totalEvents = 0;

    std::vector<std::pair<double, double>> slices; // Stores (start, end) pairs of x-slices

    int sliceEnd = startBin;
    int i = startBin;

    while( sliceEnd <= endBin ){

    // for( int i = 1; i < numBinsX; i++ ){

        while( totalEvents <= targetEvents && sliceEnd < numBinsX){
            hist->GetXaxis()->SetRange(i, sliceEnd);
            totalEvents = hist->ProjectionY()->GetEntries();
            if( print ){
                cout << "--------------" << endl;
                cout << "total events: " << totalEvents << endl;
                cout << "bins: " << i << " -> " << sliceEnd << endl;
            }

            if( totalEvents > targetEvents ){
                if( print ){
                    cout << "breaking for: " << sliceEnd << ", " << totalEvents << endl;
                }
                break;
            }   
            sliceEnd++;  
        }
        if( print ){
            cout << "_----------------------------_" << endl;
            cout << "found bin range: " << i <<  " -> " << sliceEnd  << ", Total events: " << totalEvents << endl;
            cout << "_----------------------------_" << endl;    
        }

        slices.push_back(std::make_pair(i, sliceEnd));

        totalEvents = 0;

        i = sliceEnd;
        sliceEnd = i + 1;

    }

    hist->GetXaxis()->SetRange(1, numBinsX);
    return slices;

}

#endif