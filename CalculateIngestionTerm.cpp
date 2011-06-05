#include "mex.h"
#include "CommonFunctions.hpp"
#include<cmath>
#include<cstdlib>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    unsigned no_cells = mxGetM(prhs[0]);   
    unsigned no_elements = mxGetM(prhs[2]);
    unsigned no_nodes = mxGetM(prhs[3]);
    
    double* cell_ingestion_rates = mxGetPr(prhs[1]);
    double* element_nodes = mxGetPr(prhs[2]);
    double* node_concentrations = mxGetPr(prhs[3]);
    double* previous_node_positions = mxGetPr(prhs[4]);
    double* nodes_index_in_real_nodes = mxGetPr(prhs[5]);
    double no_real_nodes = mxGetScalar(prhs[6]);

    plhs[0] = mxCreateDoubleMatrix(no_real_nodes, 1, mxREAL);
    double* ingestion_term = mxGetPr(plhs[0]);
    
	// initialise ingestion_term to zero (this may not be necessary but seems like good 
	// practise
    for(unsigned temp_FEM_node=0;temp_FEM_node<no_real_nodes;temp_FEM_node++){
        ingestion_term[temp_FEM_node] = 0;
    }
	
	// we are going to loop over every cell, then loop over the elements in that cell, and 
	// finally over the basis functions that are non-zero in each element. for each element 
	// there are three non-zero basis functions, and we will edit the corresponding
	// entries in the ingestion_term vector.
    for(unsigned current_cell_ci=0; current_cell_ci<no_cells; current_cell_ci++) {
        
//         printf("hello1 \n");
        
        double current_cell_ingestion_rate = cell_ingestion_rates[current_cell_ci];
        
        mxArray* mx_cell_elements = mxGetCell(prhs[0], current_cell_ci);
        unsigned no_cell_elements = mxGetN(mx_cell_elements);
        double* cell_elements_mi = mxGetPr(mx_cell_elements);
        
        for(unsigned current_element_local=0;current_element_local<no_cell_elements; current_element_local++){
        
//             printf("hello2 \n");
            
            unsigned current_element_mi = cell_elements_mi[current_element_local];
            unsigned current_element_ci = current_element_mi-1;
            
			// we need to store these for each for use when calculating the ingestion term.
            unsigned current_element_nodes_in_real_nodes_ci[3];
            double current_element_node_concentrations[3];
            double current_element_node_positions[3][2],jacobian_mapping[2][2];
            
			// we loop over the three basis functions that are non-zero in the current element. these 
			// are the three basis functions whose peaks are at the three nodes of the element. this is 
			// the same as looping over the three nodes of the element, it is important to remember 
			// that we are looping over basis functions
            for(unsigned current_basis_function_local=0;current_basis_function_local<3;current_basis_function_local++){
                               
                unsigned current_basis_function_mi = (unsigned)element_nodes[current_element_ci+current_basis_function_local*no_elements];
//                 printf("%i \n",current_cell_element_node_mi);
                
                unsigned current_basis_function_ci = current_basis_function_mi-1;
//                 printf("%i \n",current_cell_element_node_ci);
                
                current_element_node_concentrations[current_basis_function_local] = node_concentrations[current_basis_function_ci];
//                 printf("%f \n",current_cell_element_node_concentrations[current_basis_function_local]);
                
				// this could equally be called current_basis_function_in_real_basis_functions
                unsigned current_node_in_real_nodes_mi = nodes_index_in_real_nodes[current_basis_function_ci];
//                 printf("%i \n", current_cell_element_nodes_in_real_nodes_mi);
                
                current_element_nodes_in_real_nodes_ci[current_basis_function_local] = current_node_in_real_nodes_mi-1;
//                 printf("%i \n", current_cell_element_nodes_in_real_nodes_ci);

                for(unsigned dim=0;dim<2;dim++){
                                       
                    current_element_node_positions[current_basis_function_local][dim] = 
                            previous_node_positions[current_basis_function_ci+dim*no_nodes];
                    
                }
                
//                 printf("%f %f \n", current_cell_element_node_positions[i][0], current_cell_element_node_positions[i][1]);
            
            }
                    
            jacobian_mapping[0][0] = current_element_node_positions[2][0]-
                    current_element_node_positions[0][0];
            jacobian_mapping[0][1] = current_element_node_positions[1][0]-
                    current_element_node_positions[0][0];
            jacobian_mapping[1][0] = current_element_node_positions[2][1]-
                    current_element_node_positions[0][1];
            jacobian_mapping[1][1] = current_element_node_positions[1][1]-
                    current_element_node_positions[0][1];
            
            double det_jac = jacobian_mapping[0][0]*jacobian_mapping[1][1] -
                    jacobian_mapping[1][0]*jacobian_mapping[0][1];
            
//             printf("%f \n",det_jac);
            
			// for each basis function we edit the corresponding entry in the ingestion_term vector. for each basis function
			// we must multiply by all other basis functions, integrate, and sum. however, only three of these multiplications
			// will be non-zero, as most basis functions are zero in the current element. the only three that are non-zero are 
			// the same three that we are editing. therefore we need to multiply each function by itself and the other two before integrating.
			// when multiplying a basis function by itself on the canonical triangle, we find the integral is 1/12, whereas multiplying
			// by the other two functions the integral is 1/24.
            for(unsigned current_basis_function_local=0;current_basis_function_local<3;current_basis_function_local++){
                
                ingestion_term[current_element_nodes_in_real_nodes_ci[current_basis_function_local]] +=
                        fabs(det_jac)*current_cell_ingestion_rate*1.0/24.0*(
                        current_element_node_concentrations[0]+current_element_node_concentrations[1]+
                        current_element_node_concentrations[2]+current_element_node_concentrations[current_basis_function_local]);
                
//                 printf("%f ",ingestion_term[current_cell_element_nodes_in_real_nodes_ci[current_basis_function_local]]);
                
            }
            
//             printf("\n");
        }
    }
}

        