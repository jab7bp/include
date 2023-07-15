#ifndef calc_functions
#define calc_functions

double calc_mean(double *array){
	int n = (sizeof(array)/sizeof(array[0]));
	double sum = 0.0, mean;

	for(int i = 0; i < n; i++){
		sum += array[i];
	}
	mean = sum/n;
	return mean;
}

double calc_StdDev(double *array){
	int n = sizeof(array)/sizeof(array[0]);
	cout << "n: " << n << endl;

	return 0;
}

template<typename T>
double VectorMean(std::vector<T> const& v){
	if(v.empty()){
		return 0;
	}
	return std::accumulate(v.begin(), v.end(), 0.0)/v.size();
}

double calc_gaus_error(double par_0, double par_0_err, double par_2, double par_2_err, double bin_width ){
	double gaus_error = 0.0;

	gaus_error = ( sqrt( TMath::Pi()*( pow(par_2, 2)*pow(par_0_err, 2) + pow(par_0, 2)*pow(par_2_err, 2) ) ) )/bin_width;

	return gaus_error;
}

void list_files(TString directory, vector<TString> &filenames, TString ext){
	const char *dirname = directory.Data();

	TSystemDirectory dir(dirname, dirname);
	TList *files = dir.GetListOfFiles();

	if(files){
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while( (file = (TSystemFile*)next() ) ){
			fname = file->GetName();
			if( !file->IsDirectory() && fname.BeginsWith( ext.Data() )) {
				filenames.push_back( Form("%s/%s", directory.Data(), fname.Data()) );
				// cout << fname.Data() << endl;
			}
		}
	}
}

#endif