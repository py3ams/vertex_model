#include "mex.h"
#include<cmath>
#include<cstdlib>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
	unsigned no_cells = mxGetM(prhs[0]);   
	// N.B. This input argument must contain the nodes of the FEM_elements - not just the
	// indices. It is an Nx3 array!
    double* element_nodes = mxGetPr(prhs[1]);
	double* node_concentrations = mxGetPr(prhs[2]);
    double* node_positions = mxGetPr(prhs[3]);
    
	unsigned no_elements = mxGetM(prhs[1]);
    unsigned no_nodes = mxGetM(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix(no_cells, 1, mxREAL);
    double* cell_concentrations = mxGetPr(plhs[0]);
	
	// initialise cell_concentrations to 0.
	for(unsigned temp_cell=0;temp_cell<no_cells;temp_cell++){
		cell_concentrations[temp_cell] = 0;
	}
    
	for(unsigned current_cell_ci=0; current_cell_ci<no_cells; current_cell_ci++) {
	
		mxArray* mx_cell_elements = mxGetCell(prhs[0], current_cell_ci);
		unsigned no_cell_elements = mxGetN(mx_cell_elements);
		double* cell_elements_mi = mxGetPr(mx_cell_elements);
		
		for(unsigned current_element_local=0;current_element_local<no_cell_elements; current_element_local++){
			
			unsigned current_element_mi = cell_elements_mi[current_element_local];
			unsigned current_element_ci = current_element_mi-1;
        
			double current_element_node_concentrations[3];
			double current_element_node_positions[3][2];
			
			// loop over the three basis functions that are non-zero in the current element
			for(unsigned current_basis_function_local=0;current_basis_function_local<3;current_basis_function_local++){
				
				unsigned current_basis_function_mi = (unsigned)element_nodes[current_element_ci+current_basis_function_local*no_elements];
				unsigned current_basis_function_ci = current_basis_function_mi-1;
				
				current_element_node_concentrations[current_basis_function_local] = node_concentrations[current_basis_function_ci];
				
//             printf("%f ",Dpp_at_nodes_current_element[current_basis_function_local]);
				
				for(unsigned dim=0;dim<2;dim++){
					
					current_element_node_positions[current_basis_function_local][dim] =
							node_positions[current_basis_function_ci+dim*no_nodes];
					
				}
			}
			
//         printf("\n");
			
			double jacobian_mapping[2][2];
			
			jacobian_mapping[0][0] = current_element_node_positions[2][0]-current_element_node_positions[0][0];
			jacobian_mapping[0][1] = current_element_node_positions[1][0]-current_element_node_positions[0][0];
			jacobian_mapping[1][0] = current_element_node_positions[2][1]-current_element_node_positions[0][1];
			jacobian_mapping[1][1] = current_element_node_positions[1][1]-current_element_node_positions[0][1];
			
			double det_jac = jacobian_mapping[0][0]*jacobian_mapping[1][1] - jacobian_mapping[1][0]*jacobian_mapping[0][1];
			
			cell_concentrations[current_cell_ci] += fabs(det_jac)*1.0/6.0*(current_element_node_concentrations[0] +
					current_element_node_concentrations[1] + current_element_node_concentrations[2]);
			
			
		}
	}
}
