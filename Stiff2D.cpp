#include "mex.h"
#include "CommonFunctions.hpp"
#include<cmath>
#include<cstdlib>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	if(nrhs != 4) {
		mexErrMsgTxt("4 input arguments required");
	}
	
	double delta_t = mxGetScalar(prhs[0]);
	double* FEM_elements = mxGetPr(prhs[1]);
	double* node_positions = mxGetPr(prhs[2]);
	double* previous_node_positions = mxGetPr(prhs[3]);
	
	unsigned no_elements = mxGetM(prhs[1]);
	unsigned array_sizes = mxGetM(prhs[2]);
	
	plhs[0] = mxCreateDoubleMatrix(9*no_elements, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(9*no_elements, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(9*no_elements, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(9*no_elements, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(9*no_elements, 1, mxREAL);
	plhs[5] = mxCreateDoubleMatrix(9*no_elements, 1, mxREAL);
	plhs[6] = mxCreateDoubleMatrix(no_elements, 1, mxREAL);
	
	double* row_index = mxGetPr(plhs[0]);
	double* column_index = mxGetPr(plhs[1]);
	double* AV = mxGetPr(plhs[2]);
	double* MV = mxGetPr(plhs[3]);
	double* WV = mxGetPr(plhs[4]);
	double* MVconst = mxGetPr(plhs[5]);
	double* triangle_quality = mxGetPr(plhs[6]);
	
	unsigned counter = 0;
	
	for(unsigned current_element_ci=0;current_element_ci<no_elements;current_element_ci++){
		
		double max_edge_length = 0;
		double min_edge_length = 10;
		
		unsigned nodes_current_element_mi[3];
		
		double node_positions_current_element[3][2];
		double previous_node_positions_current_element[3][2];
		double node_velocities[3][2];
		double jacobian_mesh_velocity[2][2];
		double jacobian_mapping[2][2];
		double jacobian_inverse_mapping_transposed[2][2];
		
//         double node_position_1[2], node_position_2[2], node_position_3[2];
//         double previous_node_position_1[2], previous_node_position_2[2], previous_node_position_3[2];
		
		for(unsigned current_node_local=0;current_node_local<3;current_node_local++){
			
			double current_node_position[2], clockwise_node_position[2];
			
			unsigned current_node_global_mi = FEM_elements[current_element_ci+current_node_local*no_elements];
			unsigned current_node_global_ci = current_node_global_mi-1;
			
			unsigned clockwise_node_local = (current_node_local+1)%3;
			unsigned clockwise_node_global_mi = FEM_elements[current_element_ci+clockwise_node_local*no_elements];
			unsigned clockwise_node_global_ci = clockwise_node_global_mi-1;
			
			nodes_current_element_mi[current_node_local] = current_node_global_mi;
			
			for(unsigned k=0;k<2;k++){
				
				current_node_position[k] = node_positions[current_node_global_ci+k*array_sizes];
				clockwise_node_position[k] = node_positions[clockwise_node_global_ci+k*array_sizes];
				
				node_positions_current_element[current_node_local][k] =
						  node_positions[current_node_global_ci+k*array_sizes];
				
				previous_node_positions_current_element[current_node_local][k] =
						  previous_node_positions[current_node_global_ci+k*array_sizes];
				
				node_velocities[current_node_local][k] =
						  (node_positions_current_element[current_node_local][k] -
						  previous_node_positions_current_element[current_node_local][k])/delta_t;
				
			}
			
			double edge_length = findStraightLineDistanceBetweenTwoNodes(current_node_position, clockwise_node_position);
			
			if(edge_length>max_edge_length){
				max_edge_length = edge_length;
			}
			if(edge_length<min_edge_length){
				min_edge_length = edge_length;
			}
			
		}
		
		triangle_quality[current_element_ci] = max_edge_length/min_edge_length;
		
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
		
		jacobian_inverse_mapping_transposed[0][0] = jacobian_mapping[1][1]/det_jac;
		jacobian_inverse_mapping_transposed[0][1] = -jacobian_mapping[1][0]/det_jac;
		jacobian_inverse_mapping_transposed[1][0] = -jacobian_mapping[0][1]/det_jac;
		jacobian_inverse_mapping_transposed[1][1] = jacobian_mapping[0][0]/det_jac;
		
		double grad_phi_E[3][2] = {{-1, -1}, {0, 1}, {1, 0}};
		double grad_phi_e[3][2] = {0};
		
		for(unsigned i=0;i<3;i++){
			for(unsigned j=0;j<2;j++){
				for(unsigned k=0;k<2;k++){
					grad_phi_e[i][j] += jacobian_inverse_mapping_transposed[j][k]*grad_phi_E[i][k];
				}
			}
		}	
		
		for(unsigned i=0;i<3;i++){
			for(unsigned j=0;j<3;j++){
				
                // this is important - i is the column and j is the row!!
                row_index[counter] = nodes_current_element_mi[j];
				column_index[counter] = nodes_current_element_mi[i];
								
				if(i==j){
					MV[counter] = fabs(det_jac)*1.0/12.0;
				}
				else{
					MV[counter] = fabs(det_jac)*1.0/24.0;
				}
				
				MVconst[counter] = fabs(det_jac)*1.0/2.0;
				
				double grad_phi_i_dot_grad_phi_j = 0;
				
				for(unsigned k=0;k<2;k++){
					grad_phi_i_dot_grad_phi_j += grad_phi_e[i][k]*grad_phi_e[j][k];
				}
				AV[counter] = fabs(det_jac)*grad_phi_i_dot_grad_phi_j*1.0/2.0;
				
				WV[counter] = 0;
				for(unsigned k=0;k<2;k++){
					WV[counter] += fabs(det_jac)*grad_phi_e[j][k]*1.0/24.0*(node_velocities[0][k] +
							  node_velocities[1][k] + node_velocities[2][k] + node_velocities[i][k]);
				}
				counter++;
			}
		}
	}
}
