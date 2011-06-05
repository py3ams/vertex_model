#include<cmath>
#include<vector>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
				
	if(mxGetM(prhs[nrhs-1])>1 || mxGetN(prhs[nrhs-1])>1){
		mexErrMsgTxt("Final argument should be the required number of stats");
	}
	
	unsigned no_variables = (unsigned)(nrhs-1);
	unsigned required_no_stats = (unsigned)mxGetScalar(prhs[no_variables]);	
	
	if(required_no_stats>5){
		mexErrMsgTxt("Only 5 stats currently available");
	}
	
	plhs[0] = mxCreateDoubleMatrix(no_variables, required_no_stats, mxREAL);
	
	double* stats = mxGetPr(plhs[0]);
	double totals[no_variables];
	
	for(unsigned current_variable_index=0;current_variable_index<no_variables;current_variable_index++){
		
		double* values = mxGetPr(prhs[current_variable_index]);
		unsigned no_values = mxGetM(prhs[current_variable_index]);
		
		double max_current_variable = values[0];
		double min_current_variable = values[0];
		
		double total_current_variable = 0;
		
		for(unsigned current_value_index=0;current_value_index<no_values;current_value_index++){
			
			double current_value = values[current_value_index];
			total_current_variable += current_value;
			max_current_variable = std::max(max_current_variable, current_value);
			min_current_variable = std::min(min_current_variable, current_value);
			
		}
		
		totals[current_variable_index] = total_current_variable;
		stats[current_variable_index] = total_current_variable/(double)no_values;
		
		if(required_no_stats>1){
			
			stats[current_variable_index+no_variables] = max_current_variable;
			
			if(required_no_stats>2){
				
				stats[current_variable_index+2*no_variables] = min_current_variable;
			}
		}
	}
	
	if(required_no_stats>3){
		
		for(unsigned current_variable_index=0;current_variable_index<no_variables;current_variable_index++){
			
			double* values = mxGetPr(prhs[current_variable_index]);
			unsigned no_values = mxGetM(prhs[current_variable_index]);
			
			double mean_current_variable = stats[current_variable_index];
			double std_sum_current_variable = 0;
			
			for(unsigned current_value_index=0;current_value_index<no_values;current_value_index++){
				
				std_sum_current_variable +=
						  pow(values[current_value_index]-mean_current_variable, 2);
			}
			stats[current_variable_index+3*no_variables] =
					  sqrt(std_sum_current_variable/((double)no_values-1));
			
			if(required_no_stats>4){
				stats[current_variable_index+4*no_variables] = totals[current_variable_index];
			}
			
		}
	}
}


