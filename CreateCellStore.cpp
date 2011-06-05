#include "mex.h"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	if(nrhs != 2) {
		mexErrMsgTxt("2 input arguments required");
	}
		
	unsigned no_cells = mxGetM(prhs[0]);
	unsigned array_sizes = mxGetScalar(prhs[1]);
	
	plhs[0] = mxCreateDoubleMatrix(array_sizes, 10, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(array_sizes, 1, mxREAL);
	
	double* cell_store = mxGetPr(plhs[0]);
	double* cells_per_node = mxGetPr(plhs[1]);
	
	for(unsigned current_cell=0;current_cell<no_cells;current_cell++) {
		
		mxArray* mx_cell_nodes = mxGetCell(prhs[0], current_cell);
		
		double* cell_nodes = mxGetPr(mx_cell_nodes);
		unsigned no_cell_nodes = mxGetN(mx_cell_nodes);
		
		for(unsigned current_node_local=0;current_node_local<no_cell_nodes;current_node_local++) {
			
			unsigned current_node_global = (unsigned)cell_nodes[current_node_local];
			cell_store[current_node_global-1+(unsigned)array_sizes*
					  (unsigned)cells_per_node[current_node_global-1]] = current_cell+1;
			
			cells_per_node[current_node_global-1]++;
			
		}
	}
}
