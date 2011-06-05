#include "mex.h"
#include "CommonFunctions.hpp"
#include<cmath>
#include<cstdlib>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    unsigned no_cells = mxGetM(prhs[0]);   
    unsigned no_FEM_elements = mxGetM(prhs[3]);
    unsigned no_FEM_nodes = mxGetM(prhs[4]);
    
    double* cell_source_rates = mxGetPr(prhs[1]);
    double* cell_areas = mxGetPr(prhs[2]);
	double* element_nodes = mxGetPr(prhs[3]);
    double* previous_node_positions = mxGetPr(prhs[4]);
    double* nodes_index_in_real_nodes = mxGetPr(prhs[5]);
    double no_real_nodes = mxGetScalar(prhs[6]);

    plhs[0] = mxCreateDoubleMatrix(no_real_nodes, 1, mxREAL);
    double* source_term = mxGetPr(plhs[0]);
    
	// initialise source_term to zero (this may not be necessary but seems like good 
	// practise
    for(unsigned temp_node=0;temp_node<no_real_nodes;temp_node++){
        source_term[temp_node] = 0;
    }
	
	// we are going to loop over every cell, then loop over the elements in that cell, and 
	// finally over the basis functions that are non-zero in each element. for each element 
	// there are three non-zero basis functions, and we will edit the corresponding
	// entries in the source_term vector.
    for(unsigned current_cell_ci=0; current_cell_ci<no_cells; current_cell_ci++) {
        
//         printf("hello1 \n");
        
        double current_cell_source_rate = cell_source_rates[current_cell_ci];
        double current_cell_area = cell_areas[current_cell_ci];
		
        mxArray* mx_cell_elements = mxGetCell(prhs[0], current_cell_ci);
        unsigned no_cell_elements = mxGetN(mx_cell_elements);
        double* cell_elements_mi = mxGetPr(mx_cell_elements);
        
        for(unsigned current_element_local=0;current_element_local<no_cell_elements; current_element_local++){
        
//             printf("hello2 \n");
            
            unsigned current_element_mi = cell_elements_mi[current_element_local];
            unsigned current_element_ci = current_element_mi-1;
            
			// we need to store these for each for use when calculating the source term.
            unsigned current_element_nodes_in_real_nodes_ci[3];
            double current_element_node_concentrations[3];
            double current_element_node_positions[3][2],jacobian_mapping[2][2];
            
			// we loop over the three basis functions that are non-zero in the current element. these 
			// are the three basis functions whose peaks are at the three nodes of the element. this is 
			// the same as looping over the three nodes of the element, it is important to remember 
			// that we are looping over basis functions
            for(unsigned current_basis_function_local=0;current_basis_function_local<3;current_basis_function_local++){
                               
                unsigned current_basis_function_mi = (unsigned)element_nodes[current_element_ci+current_basis_function_local*no_FEM_elements];
//                 printf("%i \n",current_cell_element_node_mi);
                
                unsigned current_basis_function_ci = current_basis_function_mi-1;
//                 printf("%i \n",current_cell_element_node_ci);
                
				// this could equally be called current_basis_function_in_real_basis_functions
                unsigned current_node_in_real_nodes_mi = nodes_index_in_real_nodes[current_basis_function_ci];
//                 printf("%i \n", current_cell_element_nodes_in_real_nodes_mi);
                
                current_element_nodes_in_real_nodes_ci[current_basis_function_local] = current_node_in_real_nodes_mi-1;
//                 printf("%i \n", current_cell_element_nodes_in_real_nodes_ci);

                for(unsigned dim=0;dim<2;dim++){
                                       
                    current_element_node_positions[current_basis_function_local][dim] = 
                            previous_node_positions[current_basis_function_ci+dim*no_FEM_nodes];
                    
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
            
			// for each basis function we edit the corresponding entry in the source_term vector. 
			// 1/6 is the integral of a basis function over the canonical triangle.
            for(unsigned current_basis_function_local=0;current_basis_function_local<3;current_basis_function_local++){
                
                source_term[current_element_nodes_in_real_nodes_ci[current_basis_function_local]] +=
                        fabs(det_jac)*1.0/6.0*current_cell_source_rate/current_cell_area;
                
//                 printf("%f ",source_term[current_cell_element_nodes_in_real_nodes_ci[current_basis_function_local]]);
                
            }
            
//             printf("\n");
        }
    }
}

        