#include "mex.h"
#include<cmath>
#include<cstdlib>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double* cell_centroid = mxGetPr(prhs[0]);
    double* cell_nodes = mxGetPr(prhs[1]);
    double* Dpp = mxGetPr(prhs[2]);
    double* node_positions = mxGetPr(prhs[3]);
    double* total_Dpp_in_cell = mxGetPr(prhs[4]);
    
    unsigned no_cell_nodes = mxGetN(prhs[1]);
    unsigned array_sizes = mxGetM(prhs[3]);
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double* centroid_Dpp = mxGetPr(plhs[0]);
    
    double known_Dpp_total = 0;
    double unknown_Dpp_total = 0;
    
    // loop over all cell_nodes (in order to loop over elements)
    for(unsigned current_node_local=0;current_node_local<no_cell_nodes;current_node_local++){
        
        unsigned current_node_global_mi = (unsigned)cell_nodes[current_node_local];
        unsigned clockwise_node_local = (current_node_local+1)%no_cell_nodes;
        unsigned clockwise_node_global_mi = (unsigned)cell_nodes[clockwise_node_local];
        
        unsigned current_node_global_ci = current_node_global_mi-1;
        unsigned clockwise_node_global_ci = clockwise_node_global_mi-1;
        
        double Dpp_at_cell_nodes_current_element[2];
        double node_positions_current_element[3][2];
        
        Dpp_at_cell_nodes_current_element[0] = Dpp[current_node_global_ci];
        Dpp_at_cell_nodes_current_element[1] = Dpp[clockwise_node_global_ci];
        
        for(unsigned dim=0;dim<2;dim++){
            
            node_positions_current_element[0][dim] =
                    node_positions[current_node_global_ci+dim*array_sizes];
            
            node_positions_current_element[1][dim] =
                    node_positions[clockwise_node_global_ci+dim*array_sizes];
            
            node_positions_current_element[2][dim] =
                    cell_centroid[dim];
            
        }
        
        double jacobian_mapping[2][2];
        
        jacobian_mapping[0][0] = node_positions_current_element[2][0]-node_positions_current_element[0][0];
        jacobian_mapping[0][1] = node_positions_current_element[1][0]-node_positions_current_element[0][0];
        jacobian_mapping[1][0] = node_positions_current_element[2][1]-node_positions_current_element[0][1];
        jacobian_mapping[1][1] = node_positions_current_element[1][1]-node_positions_current_element[0][1];
        
        double det_jac = jacobian_mapping[0][0]*jacobian_mapping[1][1] - jacobian_mapping[1][0]*jacobian_mapping[0][1];
        
        known_Dpp_total += fabs(det_jac)*1.0/6.0*(Dpp_at_cell_nodes_current_element[0] +
                Dpp_at_cell_nodes_current_element[1]);
        
        unknown_Dpp_total += fabs(det_jac)*1.0/6.0;
        
    }
    
    *centroid_Dpp = (*total_Dpp_in_cell-known_Dpp_total)/unknown_Dpp_total;
    
}
