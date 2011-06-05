#include "mex.h"
#include<cmath>
#include<cstdlib>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double* FEM_elements = mxGetPr(prhs[0]);
    double* node_positions = mxGetPr(prhs[1]);
    
    unsigned no_elements = mxGetM(prhs[0]);
    unsigned array_sizes = mxGetM(prhs[1]);
    
    plhs[0] = mxCreateDoubleMatrix(9*no_elements, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(9*no_elements, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(9*no_elements, 1, mxREAL);

    double* row_index = mxGetPr(plhs[0]);
    double* column_index = mxGetPr(plhs[1]);
    double* MV = mxGetPr(plhs[2]);

    unsigned counter = 0;
    
    for(unsigned current_element_ci=0;current_element_ci<no_elements;current_element_ci++){
        
        unsigned nodes_current_element_mi[3];
        
        double node_positions_current_element[3][2];
        double jacobian_mapping[2][2];
               
        for(unsigned current_node_local=0;current_node_local<3;current_node_local++){
            
            unsigned current_node_global_mi = FEM_elements[current_element_ci+current_node_local*no_elements];
            unsigned current_node_global_ci = current_node_global_mi-1;
            
            nodes_current_element_mi[current_node_local] = current_node_global_mi;
            
            for(unsigned k=0;k<2;k++){
                
                node_positions_current_element[current_node_local][k] =
                        node_positions[current_node_global_ci+k*array_sizes];
                
            }
        }
        
        jacobian_mapping[0][0] = node_positions_current_element[2][0]-
                node_positions_current_element[0][0];
        jacobian_mapping[0][1] = node_positions_current_element[1][0]-
                node_positions_current_element[0][0];
        jacobian_mapping[1][0] = node_positions_current_element[2][1]-
                node_positions_current_element[0][1];
        jacobian_mapping[1][1] = node_positions_current_element[1][1]-
                node_positions_current_element[0][1];
        
        double det_jac = jacobian_mapping[0][0]*jacobian_mapping[1][1] -
                jacobian_mapping[1][0]*jacobian_mapping[0][1];
               
        for(unsigned i=0;i<3;i++){
            for(unsigned j=0;j<3;j++){
                
                row_index[counter] = nodes_current_element_mi[j];
				column_index[counter] = nodes_current_element_mi[i];
                
                if(i==j){
                    MV[counter] = fabs(det_jac)*1.0/12.0;
                }
                else{
                    MV[counter] = fabs(det_jac)*1.0/24.0;
                }
                
                counter++;
            }
        }
    }
}
