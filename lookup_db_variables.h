#ifndef LOOKUP_DB_VARIABLES_H
#define LOOKUP_DB_VARIABLES_H

TString GetLine(string filename, int db_line){
  fstream myfile (filename);
  string y;
  for(int lineno = 0; getline (myfile, y) && lineno < db_line;lineno++){
    if(lineno == db_line-1){
      break;
    }
  }
  return y;

}

TString gem_dbfilename = "db_bb.gem.dat";
TString gem_dbfile, db_dir;
string db_cut_variable = "bb.gem.crosstalk_analysis";
vector<TString> split_var;

string db_thresh_variable = "bb.gem.xtalk_ratio_threshold";
vector<TString> thresh_split_var;

string db_correct_adc_flag = "bb.gem.correct_xtalk_ADCs";
vector<TString> correct_adc_flag_split_var;


int crosstalk = -1;
double thresh = -1;
int correct_adc_flag = -1;
char* db_path;

int lookup_crosstalk(int runnum){
	db_path = getenv("DB_DIR");
	std::cout << std::endl << "Runnum for DB lookup: " << runnum << std::endl;
	std::cout << "DB dir base directory in lookup_db_variables.h ---> " << db_path << std::endl << std::endl;
	
	if( runnum < 12066 ){
		std::cout << "Runnum " << runnum << " corresponds to DB_DIR: " << Form("%s", db_path) << std::endl;
		db_dir = Form("%s", db_path);
	}
	if( runnum >= 12066 && runnum < 13239 ){
		std::cout << "Runnum " << runnum << " corresponds to DB_DIR: " << Form("%s/20211117", db_path) << std::endl;
		db_dir = Form("%s/20211117", db_path);
	}
	if( runnum >= 13239 ){
		std::cout << "Runnum " << runnum << " corresponds to DB_DIR: " << Form("%s/20220106", db_path) << std::endl;
		db_dir = Form("%s/20220106", db_path);
	}

	cout << "---------------------------------" << endl;
	cout << "Looking up Crosstalk flag in DB file: " << endl;
	cout << db_dir << endl;
	cout << "---------------------------------" << endl;
	// cout << "-------------------------" << endl;
	// cout << "Root DB path: " << db_path << endl;
	// cout << "-------------------------" << endl;
	// cout << "Final DB dir: " << db_dir.Data() << endl;
	// cout << "-------------------------" << endl;

	gem_dbfile = Form("%s/%s", db_dir.Data(), gem_dbfilename.Data());
	// string db_cut_variable = "bb.gem.maxstrip_t0";

	string filename = gem_dbfile.Data();
	fstream db_file(gem_dbfile.Data());
	string x;
	bool found_cut = false;
	int db_cut_line=1;

	if (db_file.is_open()){
	while( getline(db_file,x) ){
	  if(x.find(db_cut_variable, 0) != string::npos){
	    // cout << "found bb.gem.maxstrip_t0 on line: " << db_cut_line << endl;
	    found_cut = true;
	    if(found_cut){break;}
	  }
	db_cut_line++;
	}

	}
	string cut_variable;

	TString var = GetLine(filename, db_cut_line);

	// cout <<  "var: " << var << endl;

	TObjArray *split_var = var.Tokenize(" ");
	crosstalk = atoi(((TObjString *)(split_var->At(2)))->String());

	// cout << "Crosstalk On/Off flag from DB file: " << crosstalk << endl;
	// cout << "-------------------------" << endl << endl;
	cout << "Crosstalk flag in DB file: " << crosstalk << endl;
	cout << "-------------------------" << endl << endl;
	return crosstalk;
}

TString lookup_thresh(int runnum){
	char* db_path = getenv("DB_DIR");

	if( runnum < 12066 ){
	db_dir = Form("%s", db_path);
	}
	if( runnum <= 12066 && runnum < 13239 ){
	db_dir = Form("%s/20211117", db_path);
	}
	if( runnum >= 13239 ){
	db_dir = Form("%s/20220106", db_path);
	}

	cout << "---------------------------------" << endl;
	cout << "Looking up Xtalk Thresh in DB file: " << endl;
	cout << db_dir << endl;
	cout << "---------------------------------" << endl;
	// cout << "-------------------------" << endl;
	// cout << "Root DB path: " << db_path << endl;
	// cout << "-------------------------" << endl;
	// cout << "Final DB dir: " << db_dir.Data() << endl;
	// cout << "-------------------------" << endl;

	gem_dbfile = Form("%s/%s", db_dir.Data(), gem_dbfilename.Data());
	// string db_cut_variable = "bb.gem.maxstrip_t0";

	string filename = gem_dbfile.Data();
	fstream db_thresh_file(gem_dbfile.Data());
	string y;
	bool found_thresh = false;
	int db_thresh_line=1;

	if (db_thresh_file.is_open()){
	while( getline(db_thresh_file,y) ){
	  if(y.find(db_thresh_variable, 0) != string::npos){
	    // cout << "found bb.gem.maxstrip_t0 on line: " << db_cut_line << endl;
	    found_thresh = true;
	    if(found_thresh){break;}
	  }
	db_thresh_line++;
	}

	}
	string thresh_variable;

	TString thresh_var = GetLine(filename, db_thresh_line);

	// cout <<  "var: " << var << endl;

	TObjArray *thresh_split_var = thresh_var.Tokenize(" ");
	thresh = atof(((TObjString *)(thresh_split_var->At(2)))->String());

	// cout << "Crosstalk On/Off flag from DB file: " << crosstalk << endl;
	// cout << "-------------------------" << endl << endl;

	TString thresh_string;
	thresh_string = Form("%f", thresh);
	thresh_string.ReplaceAll(".", "_");
	thresh_string = thresh_string(0, 3);

	cout << "Xtalk threshold: " << thresh_string.Data() << endl;

	return thresh_string;
}

int lookup_correct_ADCs_flag(int runnum){
	char* db_path = getenv("DB_DIR");

	if( runnum < 12066 ){
	db_dir = Form("%s", db_path);
	}
	if( runnum <= 12066 && runnum < 13239 ){
	db_dir = Form("%s/20211117", db_path);
	}
	if( runnum >= 13239 ){
	db_dir = Form("%s/20220106", db_path);
	}

	cout << "---------------------------------" << endl;
	cout << "Looking up Correct ADC Flag in DB file: " << endl;
	cout << db_dir << endl;
	cout << "---------------------------------" << endl;
	// cout << "-------------------------" << endl;
	// cout << "Root DB path: " << db_path << endl;
	// cout << "-------------------------" << endl;
	// cout << "Final DB dir: " << db_dir.Data() << endl;
	// cout << "-------------------------" << endl;

	gem_dbfile = Form("%s/%s", db_dir.Data(), gem_dbfilename.Data());
	// string db_cut_variable = "bb.gem.maxstrip_t0";

	string filename = gem_dbfile.Data();
	fstream db_correct_adc_flag_file(gem_dbfile.Data());
	string y;
	bool found_correct_adc_flag = false;
	int db_correct_adc_flag_line=1;

	if (db_correct_adc_flag_file.is_open()){
	while( getline(db_correct_adc_flag_file,y) ){
	  if(y.find(db_correct_adc_flag, 0) != string::npos){
	    // cout << "found bb.gem.maxstrip_t0 on line: " << db_cut_line << endl;
	    found_correct_adc_flag = true;
	    if(found_correct_adc_flag){break;}
	  }
	db_correct_adc_flag_line++;
	}

	}
	string correct_adc_flag_variable;

	TString correct_adc_flag_var = GetLine(filename, db_correct_adc_flag_line);

	// cout <<  "var: " << var << endl;

	TObjArray *correct_adc_flag_split_var = correct_adc_flag_var.Tokenize(" ");
	correct_adc_flag = atoi(((TObjString *)(correct_adc_flag_split_var->At(2)))->String());

	// cout << "Crosstalk On/Off flag from DB file: " << crosstalk << endl;
	// cout << "-------------------------" << endl << endl;

	TString correct_adc_flag_string;
	correct_adc_flag_string = Form("%d", correct_adc_flag);
	correct_adc_flag_string.ReplaceAll(".", "_");
	correct_adc_flag_string = correct_adc_flag_string(0, 3);

	cout << "Xtalk threshold: " << correct_adc_flag_string.Data() << endl;

	return correct_adc_flag;
}

#endif