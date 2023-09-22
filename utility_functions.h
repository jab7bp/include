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


#endif