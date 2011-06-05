#include "mex.h"
#include<cmath>
#include<vector>
#include<cassert>
#include<cstdlib>
#include "CommonFunctions.hpp"

#define pi 3.1415926535897932384626433832795

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	if(nrhs != 5) {
		mexErrMsgTxt("Exactly 5 input arguments required");
	}
	
	//	Get array_sizes
	const mxArray* mx_array_sizes = mexGetVariablePtr("global", "array_sizes");
	if (mx_array_sizes == NULL) {
		mexErrMsgTxt("array_sizes not defined");
	}
	unsigned array_sizes = (unsigned)mxGetScalar(mx_array_sizes);
	
	//	Get iteration
	const mxArray* mx_iteration = mexGetVariablePtr("caller", "iteration");
	if (mx_iteration == NULL) {
		mexErrMsgTxt("iteration not defined");
	}
	double iteration = mxGetScalar(mx_iteration);
	
	//	Get join_probability
	const mxArray* mx_join_probability = mexGetVariablePtr("caller", "join_probability");
	if (mx_join_probability == NULL) {
		mexErrMsgTxt("join_probability not defined");
	}
	double join_probability = mxGetScalar(mx_join_probability);
	
	//	Get protection_time
	const mxArray* mx_protection_time = mexGetVariablePtr("caller", "protection_time");
	if (mx_protection_time == NULL) {
		mexErrMsgTxt("protection_time not defined");
	}
	double protection_time = mxGetScalar(mx_protection_time);
	
	//	Get threshold_join_nodes
	const mxArray* mx_threshold_join_nodes = mexGetVariablePtr("caller", "threshold_join_nodes");
	if (mx_threshold_join_nodes == NULL) {
		mexErrMsgTxt("threshold_join_nodes not defined");
	}
	double threshold_join_nodes = mxGetScalar(mx_threshold_join_nodes);
	
	//	Get time
	const mxArray* mx_time = mexGetVariablePtr("global", "time");
	if (mx_time == NULL) {
		mexErrMsgTxt("time not defined");
	}
	double time = mxGetScalar(mx_time);
	
	unsigned no_cells = mxGetM(prhs[1]);
	unsigned width_cell_store = mxGetN(prhs[0]);
	
	double* initial_cell_store = mxGetPr(prhs[0]);
	double* initial_cells_per_node = mxGetPr(prhs[2]);
	double* initial_node_positions = mxGetPr(prhs[3]);
	double* initial_time_nodes_created = mxGetPr(prhs[4]);
	
	plhs[0] = mxCreateDoubleMatrix(array_sizes, width_cell_store, mxREAL);
	plhs[1] = mxCreateCellMatrix(no_cells, 1);
	plhs[2] = mxCreateDoubleMatrix(array_sizes, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(array_sizes, 2, mxREAL);
	plhs[5] = mxCreateDoubleMatrix(array_sizes, 1, mxREAL);
	
	double* final_cell_store = mxGetPr(plhs[0]);
	double* final_cells_per_node = mxGetPr(plhs[2]);
	double* no_join_nodes_this_iteration = mxGetPr(plhs[3]);
	double* final_node_positions = mxGetPr(plhs[4]);
	double* final_time_nodes_created = mxGetPr(plhs[5]);
	
	*no_join_nodes_this_iteration = 0;
	
	for(unsigned current_node_global_ci=0;current_node_global_ci<array_sizes;current_node_global_ci++){
		final_cells_per_node[current_node_global_ci] = initial_cells_per_node[current_node_global_ci];
		final_time_nodes_created[current_node_global_ci] = initial_time_nodes_created[current_node_global_ci];
		for(unsigned dim=0;dim<2;dim++){
			final_node_positions[current_node_global_ci+dim*array_sizes] = initial_node_positions[current_node_global_ci+array_sizes*dim];
		}
		for(unsigned i=0;i<width_cell_store;i++){
			final_cell_store[current_node_global_ci+i*array_sizes] = initial_cell_store[current_node_global_ci+i*array_sizes];
		}
	}
	
	for(unsigned current_cell=0; current_cell<no_cells; current_cell++) {
		
		mxArray* mx_initial_cell_nodes = mxGetCell(prhs[1], current_cell);
		double* initial_cell_nodes = mxGetPr(mx_initial_cell_nodes);
		unsigned no_cell_nodes = mxGetN(mx_initial_cell_nodes);
		
		mxArray* mx_final_cell_nodes = mxCreateDoubleMatrix(1, no_cell_nodes, mxREAL);
		double* final_cell_nodes = mxGetPr(mx_final_cell_nodes);
		
		for(unsigned i=0; i<no_cell_nodes;i++){
			final_cell_nodes[i] = initial_cell_nodes[i];
		}
		mxSetCell(plhs[1], current_cell, mx_final_cell_nodes);
	}
	
	for(unsigned current_cell_ci=0; current_cell_ci<no_cells; current_cell_ci++) {
		
		unsigned current_cell_mi = current_cell_ci+1;
		
// 		printf("start of main loop \n");
		
		mxArray* mx_cell_nodes = mxGetCell(plhs[1], current_cell_ci);
		double* cell_nodes = mxGetPr(mx_cell_nodes);
		
// 		for(unsigned i=0;i<6;i++){
// 			printf("%f \n",cell_nodes[i]);
// 		}
		
		unsigned no_cell_nodes = mxGetN(mx_cell_nodes);
// 		printf("%u \n",no_cell_nodes);
		
		if(no_cell_nodes > 3){
			
			for(unsigned current_node_local=0; current_node_local<no_cell_nodes; current_node_local++) {
				
				double rand_number = rand()/((double)RAND_MAX);
				//                 printf("%f \n",rand_number);
				
				if(rand_number < join_probability){
					
					unsigned current_node_global_mi = (unsigned)cell_nodes[current_node_local];
					unsigned current_node_global_ci = current_node_global_mi-1;
					unsigned clockwise_node_local;
					
// 				printf("%f \n",final_cells_per_node[current_node_global-1]);
					
					if(current_node_local==no_cell_nodes-1) {
						clockwise_node_local = 0;
					}
					else {
						clockwise_node_local = current_node_local+1;
					}
					unsigned clockwise_node_global_mi = (unsigned)cell_nodes[clockwise_node_local];
					unsigned clockwise_node_global_ci = clockwise_node_global_mi - 1;
					
					// extract the positions of the three nodes in question from the initial node positions matrix.
					// have to subtract 1 from each value as they are in matlab indices.
					double current_node_position[2], clockwise_node_position[2];
					
					for(unsigned dim=0;dim<2;dim++) {
						
						current_node_position[dim] = final_node_positions[current_node_global_ci + dim*array_sizes];
						clockwise_node_position[dim] = final_node_positions[clockwise_node_global_ci + dim*array_sizes];
					}
					
					// find current edge length (distance from current node to anti clockwise node)
					double current_edge_length = findStraightLineDistanceBetweenTwoNodes(current_node_position,
							  clockwise_node_position);
					
// 				printf("hello1 \n");
					
					if(current_edge_length < threshold_join_nodes &&
							  (time-final_time_nodes_created[current_node_global_ci])>protection_time &&
							  (time-final_time_nodes_created[clockwise_node_global_ci])>protection_time){
						
// 					printf("hello2 \n");
						
						unsigned no_cells_containing_just_node_1 = (unsigned)final_cells_per_node[current_node_global_ci]-2;
						unsigned no_cells_containing_just_node_2 = (unsigned)final_cells_per_node[clockwise_node_global_ci]-2;
						
						unsigned cells_containing_just_node_1_mi[no_cells_containing_just_node_1];
						unsigned cells_containing_just_node_2_mi[no_cells_containing_just_node_2];
						int cell_with_same_edge_mi = -1;
						
						unsigned counter_1 = 0;
						for(unsigned i=0;i<final_cells_per_node[current_node_global_ci];i++){
							
							unsigned temp_cell_1_mi = (unsigned)final_cell_store[current_node_global_ci + i*array_sizes];
							
							if(temp_cell_1_mi != current_cell_mi){
								for(unsigned j=0;j<final_cells_per_node[clockwise_node_global_ci];j++){
									
									unsigned temp_cell_2_mi = (unsigned)final_cell_store[clockwise_node_global_ci + j*array_sizes];
									
									if(temp_cell_2_mi == temp_cell_1_mi){
										cell_with_same_edge_mi = temp_cell_1_mi;
									}
								}
							}
							
							if(temp_cell_1_mi!=current_cell_mi && temp_cell_1_mi!=cell_with_same_edge_mi){
								cells_containing_just_node_1_mi[counter_1] = temp_cell_1_mi;
								counter_1++;
							}
						}
						
						
						unsigned cell_with_same_edge_ci = cell_with_same_edge_mi-1;
						
						unsigned counter_2 = 0;
						for(unsigned i=0;i<final_cells_per_node[clockwise_node_global_ci];i++){
							unsigned temp_cell_2_mi = (unsigned)final_cell_store[clockwise_node_global_ci + i*array_sizes];
							
							if(temp_cell_2_mi!=current_cell_mi && temp_cell_2_mi!=cell_with_same_edge_mi){
								cells_containing_just_node_2_mi[counter_2] = temp_cell_2_mi;
								counter_2++;
							}
						}
						
						mxArray* mx_cell_nodes_cell_with_same_edge = mxGetCell(plhs[1], cell_with_same_edge_ci);
						double* cell_nodes_cell_with_same_edge = mxGetPr(mx_cell_nodes_cell_with_same_edge);
						unsigned no_cell_nodes_cell_with_same_edge = mxGetN(mx_cell_nodes_cell_with_same_edge);
						
						if(no_cell_nodes_cell_with_same_edge > 3){
							
// 						printf("doing a join \n");
							
							(*no_join_nodes_this_iteration)++;
							
							unsigned new_node_ci = 0;
							while(final_cells_per_node[new_node_ci]>0){
								new_node_ci++;
							}
							unsigned new_node_mi = new_node_ci+1;
							
// 						printf("%u ",current_node_global_mi);
// 						printf("%u ",clockwise_node_global_mi);
// 						printf("%u \n",new_node_mi);
							
// 						printf("%u ",cell_with_same_edge_mi);
// 						printf("%u \n",current_cell_mi);
							
							assert(new_node_mi <= array_sizes);
							
							final_time_nodes_created[new_node_ci] = time;
							
							for(unsigned dim=0;dim<2;dim++){
								final_node_positions[new_node_ci+dim*array_sizes] =
										  (current_node_position[dim]+clockwise_node_position[dim])/2;
								
							}
							
							mxArray* mx_cell_nodes_edited = mxCreateDoubleMatrix(1, no_cell_nodes-1, mxREAL);
							double* cell_nodes_edited = mxGetPr(mx_cell_nodes_edited);
							
							int temp_counter = -1;
							for(unsigned temp_current_node_local=0;temp_current_node_local<no_cell_nodes;temp_current_node_local++){
								
								unsigned temp_current_node_global_mi = (unsigned)cell_nodes[temp_current_node_local];
								
								if(temp_current_node_global_mi==current_node_global_mi){}
								else if(temp_current_node_global_mi==clockwise_node_global_mi){
									temp_counter++;
									cell_nodes_edited[temp_counter]=new_node_mi;
								}
								else{
									temp_counter++;
									cell_nodes_edited[temp_counter]=temp_current_node_global_mi;
								}
							}
							mxSetCell(plhs[1], current_cell_ci, mx_cell_nodes_edited);
							
							mxArray* mx_cell_nodes_cell_with_same_edge_edited = mxCreateDoubleMatrix(1, no_cell_nodes_cell_with_same_edge-1, mxREAL);
							double* cell_nodes_cell_with_same_edge_edited = mxGetPr(mx_cell_nodes_cell_with_same_edge_edited);
							
							temp_counter = -1;
							for(unsigned temp_current_node_local=0;temp_current_node_local<no_cell_nodes_cell_with_same_edge;temp_current_node_local++){
								
								unsigned temp_current_node_global_mi = (unsigned)cell_nodes_cell_with_same_edge[temp_current_node_local];
								
								if(temp_current_node_global_mi==current_node_global_mi){}
								else if(temp_current_node_global_mi==clockwise_node_global_mi){
									temp_counter++;
									cell_nodes_cell_with_same_edge_edited[temp_counter]=new_node_mi;
								}
								else{
									temp_counter++;
									cell_nodes_cell_with_same_edge_edited[temp_counter]=temp_current_node_global_mi;
								}
							}
							mxSetCell(plhs[1], cell_with_same_edge_ci, mx_cell_nodes_cell_with_same_edge_edited);
							
							for(unsigned temp_cell_local=0;temp_cell_local<no_cells_containing_just_node_1;temp_cell_local++){
								
								unsigned temp_cell_global_mi = (unsigned)cells_containing_just_node_1_mi[temp_cell_local];
								unsigned temp_cell_global_ci = temp_cell_global_mi-1;
								
								mxArray* mx_temp_cell_nodes = mxGetCell(plhs[1], temp_cell_global_ci);
								double* temp_cell_nodes = mxGetPr(mx_temp_cell_nodes);
								unsigned no_temp_cell_nodes = mxGetN(mx_temp_cell_nodes);
								
								mxArray* mx_temp_cell_nodes_edited = mxCreateDoubleMatrix(1, no_temp_cell_nodes, mxREAL);
								double* temp_cell_nodes_edited = mxGetPr(mx_temp_cell_nodes_edited);
								
								for(unsigned temp_cell_node_local=0;temp_cell_node_local<no_temp_cell_nodes;temp_cell_node_local++){
									
									unsigned temp_cell_node_global_mi = (unsigned)temp_cell_nodes[temp_cell_node_local];
									
									if(temp_cell_node_global_mi == current_node_global_mi){
										temp_cell_nodes_edited[temp_cell_node_local] = new_node_mi;
									}
									else{
										temp_cell_nodes_edited[temp_cell_node_local] = temp_cell_nodes[temp_cell_node_local];
									}
								}
								
								mxSetCell(plhs[1], temp_cell_global_ci, mx_temp_cell_nodes_edited);
							}
							
							for(unsigned temp_cell_local=0;temp_cell_local<no_cells_containing_just_node_2;temp_cell_local++){
								
								unsigned temp_cell_global_mi = (unsigned)cells_containing_just_node_2_mi[temp_cell_local];
								unsigned temp_cell_global_ci = temp_cell_global_mi-1;
								
								mxArray* mx_temp_cell_nodes = mxGetCell(plhs[1], temp_cell_global_ci);
								double* temp_cell_nodes = mxGetPr(mx_temp_cell_nodes);
								unsigned no_temp_cell_nodes = mxGetN(mx_temp_cell_nodes);
								
								mxArray* mx_temp_cell_nodes_edited = mxCreateDoubleMatrix(1, no_temp_cell_nodes, mxREAL);
								double* temp_cell_nodes_edited = mxGetPr(mx_temp_cell_nodes_edited);
								
								for(unsigned temp_cell_node_local=0;temp_cell_node_local<no_temp_cell_nodes;temp_cell_node_local++){
									
									unsigned temp_cell_node_global_mi = (unsigned)temp_cell_nodes[temp_cell_node_local];
									
									if(temp_cell_node_global_mi == clockwise_node_global_mi){
										temp_cell_nodes_edited[temp_cell_node_local] = new_node_mi;
									}
									else{
										temp_cell_nodes_edited[temp_cell_node_local] = temp_cell_nodes[temp_cell_node_local];
									}
									
								}
								
								mxSetCell(plhs[1], temp_cell_global_ci, mx_temp_cell_nodes_edited);
							}
							
							final_cells_per_node[current_node_global_ci] = 0;
							final_cells_per_node[clockwise_node_global_ci] = 0;
							
							for(unsigned i=0;i<width_cell_store;i++){
								final_cell_store[current_node_global_ci+i*array_sizes] = 0;
								final_cell_store[clockwise_node_global_ci+i*array_sizes] = 0;
							}
							
							for(unsigned dim=0;dim<2;dim++){
								final_node_positions[current_node_global_ci+dim*array_sizes] = 0;
								final_node_positions[clockwise_node_global_ci+dim*array_sizes] = 0;
							}
							
							final_cells_per_node[new_node_ci] = no_cells_containing_just_node_1 +
									  no_cells_containing_just_node_2 + 2;
							
							final_cell_store[new_node_ci] = current_cell_mi;
							final_cell_store[new_node_ci+array_sizes] = cell_with_same_edge_mi;
							
							unsigned last_counter = 2;
							
							for(unsigned temp_cell_local=0;temp_cell_local<no_cells_containing_just_node_1;temp_cell_local++){
								
								unsigned temp_cell_global_mi = (unsigned)cells_containing_just_node_1_mi[temp_cell_local];
								
								final_cell_store[new_node_ci+last_counter*array_sizes] = temp_cell_global_mi;
								last_counter++;
							}
							
							for(unsigned temp_cell_local=0;temp_cell_local<no_cells_containing_just_node_2;temp_cell_local++){
								
								unsigned temp_cell_global_mi = (unsigned)cells_containing_just_node_2_mi[temp_cell_local];
								
								final_cell_store[new_node_ci+last_counter*array_sizes] = temp_cell_global_mi;
								last_counter++;
							}
							break;
						}
					}
				}
			}
		}
	}
}
