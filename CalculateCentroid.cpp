#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	// there must be 2 input arguments - cell node positions and cell area
	if (nrhs != 2) {
		mexErrMsgTxt("Exactly 2 input arguments required - cell node positions, cell area");
	}
	
	unsigned no_cell_nodes = mxGetM(prhs[0]);
		
	double* cell_node_positions = mxGetPr(prhs[0]);
	double* cell_area = mxGetPr(prhs[1]);
	
	plhs[0] = mxCreateDoubleMatrix(1, 2, mxREAL);
	double* cell_centre = mxGetPr(plhs[0]);
	
	double cx = 0, cy = 0;
	double x_i, x_iplus1, y_i, y_iplus1;
	
	for(unsigned i=0;i<no_cell_nodes;i++){
		
		x_i = cell_node_positions[i];
		x_iplus1 = cell_node_positions[(i+1)%no_cell_nodes];
		
//		printf("%f %f \n",x_i,x_iplus1);
		
		y_i = cell_node_positions[i+no_cell_nodes];
		y_iplus1 = cell_node_positions[(i+1)%no_cell_nodes+no_cell_nodes];
		
		cx = cx + (x_i + x_iplus1)*(x_i*y_iplus1 - x_iplus1*y_i);
		cy = cy + (y_i + y_iplus1)*(x_i*y_iplus1 - x_iplus1*y_i);
	}
	
	cx = cx/(6.0*(*cell_area));
	cy = cy/(6.0*(*cell_area));

	cell_centre[0] = -cx;
	cell_centre[1] = -cy;
	
}

