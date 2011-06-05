#include "mex.h"
#include<cmath>
#include<cstdlib>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double* Dpp = mxGetPr(prhs[0]);
    // N.B. This input argument must contain the nodes of the FEM_elements - not just the 
    // indices. It is an Nx3 array!
    double* FEM_elements = mxGetPr(prhs[1]);
    double* node_positions = mxGetPr(prhs[2]);
    
    unsigned no_elements = mxGetM(prhs[1]);
    unsigned no_FEM_nodes = mxGetM(prhs[2]);
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double* total_Dpp = mxGetPr(plhs[0]);
    
    *total_Dpp = 0;
	    
    // loop over all elements
    for(unsigned current_element_ci=0;current_element_ci<no_elements;current_element_ci++){
        
        double Dpp_at_nodes_current_element[3];
        double node_positions_current_element[3][2];
        
        for(unsigned current_node_local=0;current_node_local<3;current_node_local++){
            
            unsigned current_node_global_mi = (unsigned)FEM_elements[current_element_ci+current_node_local*no_elements];
            unsigned current_node_global_ci = current_node_global_mi-1;
            
            Dpp_at_nodes_current_element[current_node_local] = Dpp[current_node_global_ci];
            
//             printf("%f ",Dpp_at_nodes_current_element[current_node_local]);
            
            for(unsigned dim=0;dim<2;dim++){
                
                node_positions_current_element[current_node_local][dim] =
                        node_positions[current_node_global_ci+dim*no_FEM_nodes];
                
            }
        }
        
//         printf("\n");
        
        double jacobian_mapping[2][2];
        
        jacobian_mapping[0][0] = node_positions_current_element[2][0]-node_positions_current_element[0][0];
        jacobian_mapping[0][1] = node_positions_current_element[1][0]-node_positions_current_element[0][0];
        jacobian_mapping[1][0] = node_positions_current_element[2][1]-node_positions_current_element[0][1];
        jacobian_mapping[1][1] = node_positions_current_element[1][1]-node_positions_current_element[0][1];
        
        double det_jac = jacobian_mapping[0][0]*jacobian_mapping[1][1] - jacobian_mapping[1][0]*jacobian_mapping[0][1];
        
        *total_Dpp += fabs(det_jac)*1.0/2.0*(Dpp_at_nodes_current_element[0] +
                Dpp_at_nodes_current_element[1] + Dpp_at_nodes_current_element[2]);
		
        
    }
}
