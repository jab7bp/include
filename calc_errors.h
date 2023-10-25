#ifndef CALC_ERRORS_H
#define CALC_ERRORS_H

//Error propagation for addition and subtraction:
double CalculateErrorAdditionSubtraction( double val1, double val1_err, double val2, double val2_err){
    double error = std::sqrt(std::pow(val1_err, 2) + std::pow(val2_err, 2));
    return error;
}

//Error propagation for multiplication and division:
double CalculateErrorMultiplicationDivision( double val1, double val1_err, double val2, double val2_err, double eval ){
    double relativeError = std::sqrt(std::pow(val1_err / val1, 2) + std::pow(val2_err / val2, 2));
    double finalError = eval*relativeError;

    return finalError;
}

// Function to calculate the error in a fit function result at a given x value
double calculateFitFunctionError(TF1* fitFunction, const std::vector<double>& coefficients, const std::vector<double>& coefficientErrors, double xValue) {
    // Set the coefficients with their respective errors
    for (int i = 0; i < coefficients.size(); ++i) {
        fitFunction->SetParameter(i, coefficients[i]);
        fitFunction->SetParError(i, coefficientErrors[i]);
    }

    // Calculate the derivative of the fit function with respect to each coefficient using a finite difference method
    const double epsilon = 1e-5;  // Small finite difference step
    std::vector<double> derivatives;

    for (int i = 0; i < coefficients.size(); ++i) {
        // Calculate the forward difference
        double originalValue = fitFunction->Eval(xValue);
        fitFunction->SetParameter(i, coefficients[i] + epsilon);
        double newValue = fitFunction->Eval(xValue);
        fitFunction->SetParameter(i, coefficients[i]);

        double derivative = (newValue - originalValue) / epsilon;
        derivatives.push_back(derivative);
    }

    // Calculate the error in the result using error propagation
    double error = 0.0;
    for (int i = 0; i < coefficients.size(); ++i) {
        error += TMath::Power(derivatives[i] * coefficientErrors[i], 2);
    }

    return TMath::Sqrt(error);
}

double calculateVectorStandardDeviation(const std::vector<double>& data) {
    double mean = 0.0;
    for (double value : data) {
        mean += value;
    }
    mean /= data.size();

    double variance = 0.0;
    for (double value : data) {
        variance += (value - mean) * (value - mean);
    }
    variance /= data.size();

    return std::sqrt(variance);
}

double calculateVectorMean(const std::vector<double>& data) {
    double sum = 0.0;
    for (double value : data) {
        sum += value;
    }
    return sum / data.size();
}

double calculatePolynomialValue(const std::vector<double>& coefficients, double x) {
    double result = 0.0;
    for (int i = 0; i < coefficients.size(); ++i) {
        result += coefficients[i] * TMath::Power(x, i);
    }
    return result;
}

// Function to calculate the error using error propagation
double calculatePolynomialError(const std::vector<double>& coefficients, const std::vector<double>& coefficientErrors, double x) {
    double error = 0.0;

    // Calculate the polynomial value at x
    double polynomialValue = calculatePolynomialValue(coefficients, x);

    for (int i = 1; i < coefficients.size(); ++i) {
        // Calculate the partial derivative of the polynomial with respect to coefficient a_i
        double derivative = TMath::Power(x, i);

        // Calculate the error contribution for this coefficient
        double errorContribution = TMath::Power(derivative * coefficientErrors[i], 2);

        error += errorContribution;
    }

    // Take the square root of the sum to get the total error
    return TMath::Sqrt(error);
}

//Lets make a general error calculation function for any polynomial evaluated at x. 
double CalculateFunctionErrorAtX( vector<double> coeff, vector<double> coeff_err, double x_eval) {
    // Create a TFormula object to define your fit function (Nth order polynomial)
    //We'll get the order of the polynomial from the size of coeff and coeff_err
    TString function_str = "[0]";
    int polN = coeff.size();

    for( int i = 1; i < polN; i++ ){
        function_str = Form("%s + [%i]", function_str.Data(), i );
        for( int j = 0; j < i; j++ ){
            function_str = Form("%s*x",function_str.Data());
        }
    }
    cout << "-----------------------------------------------------" << endl;
    cout << "Evaluating the pol" << coeff.size() << " function: " << function_str.Data() << endl;
    cout << "-----------------------------------------------------" << endl;

    //Define a formula from this
    TFormula fit_function("fit_function", function_str.Data() );
    // TFormula f("fit_function", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x");
   // Set the coefficients
    for( int i = 0; i < polN; i++ ){
        fit_function.SetParameter( i, coeff[i] );
    }

    // Calculate the value of the function at the point x
    double f_x = fit_function.Eval(x_eval);

    // Calculate the partial derivatives with respect to the coefficients
    // vector<double> df_dN = {0.0, 1.0}; //First derivative of the polN
    vector<double> df_dN = {}; //First derivative of the polN (dF(x)/d(A)) is 1 but we can calculate that on next line.
    for( int power = 0; power < polN; power++ ){
        df_dN.push_back( pow(x_eval, power) );
    }

    // Calculate the error using the error propagation formula for a polynomial expansion
    double sigma_f_x = 0;

    for( int i = 0; i < polN; i++){
        sigma_f_x += pow( df_dN[i]*coeff_err[i], 2 );
    }

    // Now, sigma_f_x contains the error at the point x
    //We need the sqrt though...
    cout << "Error at " << x_eval << " = " << sqrt(sigma_f_x) << endl;

    return sqrt(sigma_f_x);
}

#endif