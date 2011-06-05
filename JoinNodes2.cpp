#include "mex.h"
#include<cmath>
#include<vector>
#include<cassert>
#include "CommonFunctions.hpp"

#define pi 3.1415926535897932384626433832795

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
// 	printf("start of file \n");
	
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
	
	//	Get protection_time
	const mxArray* mx_protection_time = mexGetVariablePtr("caller", "protection_time");
	if (mx_protection_time == NULL) {
		mexErrMsgTxt("protection_time not defined");
	}
	double protection_time = mxGetScalar(mx_protection_time);
	
	//	Get join_probability
	const mxArray* mx_join_probability = mexGetVariablePtr("caller", "join_probability");
	if (mx_join_probability == NULL) {
		mexErrMsgTxt("join_probability not defined");
	}
	double join_probability = mxGetScalar(mx_join_probability);
	
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
	
	unsigned width_cell_store = mxGetN(prhs[1]);
	unsigned no_cells = mxGetM(prhs[2]);
	unsigned no_FEM_elements = mxGetM(prhs[5]);
	unsigned length_FEM_node_positions = mxGetM(prhs[6]);
	
	double* initial_cell_store = mxGetPr(prhs[1]);
	double* initial_cells_per_node = mxGetPr(prhs[3]);
	double* initial_Dpp = mxGetPr(prhs[4]);
	double* initial_FEM_elements = mxGetPr(prhs[5]);
	double* initial_previous_FEM_node_positions = mxGetPr(prhs[6]);
	double* initial_node_positions = mxGetPr(prhs[7]);
	double* initial_time_nodes_created = mxGetPr(prhs[8]);
	
	plhs[0] = mxCreateCellMatrix(no_cells, 1);
	plhs[1] = mxCreateDoubleMatrix(array_sizes, width_cell_store, mxREAL);
	plhs[2] = mxCreateCellMatrix(no_cells, 1);
	plhs[3] = mxCreateDoubleMatrix(array_sizes, 1, mxREAL);
	plhs[4] = mxCreateDoubleMatrix(length_FEM_node_positions, 1, mxREAL);
	plhs[5] = mxCreateDoubleMatrix(no_FEM_elements, 3, mxREAL);
	plhs[6] = mxCreateDoubleMatrix(length_FEM_node_positions, 2, mxREAL);
	plhs[7] = mxCreateDoubleMatrix(array_sizes, 2, mxREAL);
	plhs[8] = mxCreateDoubleMatrix(array_sizes, 1, mxREAL);
	plhs[9] = mxCreateDoubleMatrix(1, 1, mxREAL);
	
	double* final_cell_store = mxGetPr(plhs[1]);
	double* final_cells_per_node = mxGetPr(plhs[3]);
	double* final_Dpp = mxGetPr(plhs[4]);
	double* final_FEM_elements = mxGetPr(plhs[5]);
	double* final_previous_FEM_node_positions = mxGetPr(plhs[6]);
	double* final_node_positions = mxGetPr(plhs[7]);
	double* final_time_nodes_created = mxGetPr(plhs[8]);
	double* no_join_nodes_this_iteration = mxGetPr(plhs[9]);
	
	*no_join_nodes_this_iteration = 0;
	
	// set output matrices equal to input matrices
	for(unsigned current_node=0;current_node<array_sizes;current_node++){
		final_cells_per_node[current_node] = initial_cells_per_node[current_node];
		final_time_nodes_created[current_node] = initial_time_nodes_created[current_node];
		for(unsigned dim=0;dim<2;dim++){
			final_node_positions[current_node+array_sizes*dim] = initial_node_positions[current_node+array_sizes*dim];
		}
		for(unsigned i=0;i<width_cell_store;i++){
			final_cell_store[current_node+i*array_sizes] = initial_cell_store[current_node+i*array_sizes];
		}
	}
	
	for(unsigned current_element=0;current_element<no_FEM_elements;current_element++){
		for(unsigned i=0;i<3;i++){
			final_FEM_elements[current_element+no_FEM_elements*i] =
					  initial_FEM_elements[current_element+no_FEM_elements*i];
		}
	}
	
	for(unsigned current_FEM_node=0;current_FEM_node<length_FEM_node_positions;current_FEM_node++){
		final_Dpp[current_FEM_node] = initial_Dpp[current_FEM_node];
		for(unsigned dim=0;dim<2;dim++){
			final_previous_FEM_node_positions[current_FEM_node+length_FEM_node_positions*dim] =
					  initial_previous_FEM_node_positions[current_FEM_node+length_FEM_node_positions*dim];
		}
	}
	
	// set output cells equal to input cells
	for(unsigned current_cell=0; current_cell<no_cells; current_cell++) {
		
		mxArray* mx_initial_cell_nodes = mxGetCell(prhs[2], current_cell);
		
		unsigned no_cell_nodes = mxGetN(mx_initial_cell_nodes);
		
		double* initial_cell_nodes = mxGetPr(mx_initial_cell_nodes);
		
		mxArray* mx_final_cell_nodes = mxCreateDoubleMatrix(1, no_cell_nodes, mxREAL);
		double* final_cell_nodes = mxGetPr(mx_final_cell_nodes);
		
		mxArray* mx_initial_cell_elements = mxGetCell(prhs[0], current_cell);
		double* initial_cell_elements = mxGetPr(mx_initial_cell_elements);
		
		mxArray* mx_final_cell_elements = mxCreateDoubleMatrix(1, no_cell_nodes, mxREAL);
		double* final_cell_elements = mxGetPr(mx_final_cell_elements);
		
		for(unsigned i=0; i<no_cell_nodes;i++){
			final_cell_nodes[i] = initial_cell_nodes[i];
			final_cell_elements[i] = initial_cell_elements[i];
		}
		
		mxSetCell(plhs[2], current_cell, mx_final_cell_nodes);
		mxSetCell(plhs[0], current_cell, mx_final_cell_elements);
		
	}
	
	bool join_logical;
	
	// loop over all cells
	for(unsigned current_cell_ci=0; current_cell_ci<no_cells; current_cell_ci++) {
		
		join_logical = false;
		
// 		printf("start of main loop \n");
		
		unsigned current_cell_mi = current_cell_ci+1;
		
		mxArray* mx_cell_nodes = mxGetCell(plhs[2], current_cell_ci);
		double* cell_nodes = mxGetPr(mx_cell_nodes);
		
		unsigned no_cell_nodes = mxGetN(mx_cell_nodes);
		
		// only proceed if there are more than 3 cell nodes
		if(no_cell_nodes > 3){
			
			// loop over the nodes of the current cell
			for(unsigned current_node_local=0; current_node_local<no_cell_nodes; current_node_local++) {
				
				double rand_number = rand()/((double)RAND_MAX);
				
				// only proceed with p(join_probability)
				if(rand_number < join_probability){
					
					// find current node and clockwise node in matlab and c indicies
					unsigned current_node_global_mi = (unsigned)cell_nodes[current_node_local];
					unsigned current_node_global_ci = current_node_global_mi-1;
					
					unsigned clockwise_node_local = (current_node_local+1)%no_cell_nodes;
					
					unsigned clockwise_node_global_mi = (unsigned)cell_nodes[clockwise_node_local];
					unsigned clockwise_node_global_ci = clockwise_node_global_mi - 1;
					
					// extract the positions of the current node and clockwise node from the node positions matrix
					double current_node_position[2], clockwise_node_position[2];
					
					for(unsigned dim=0;dim<2;dim++) {
						
						current_node_position[dim] = final_node_positions[current_node_global_ci + dim*array_sizes];
						clockwise_node_position[dim] = final_node_positions[clockwise_node_global_ci + dim*array_sizes];
					}
					
					// find the current edge length
					double current_edge_length = findStraightLineDistanceBetweenTwoNodes(current_node_position,
							  clockwise_node_position);
					
					// only proceed if edge length is less than threshold, and nodes have not been created very
					// recently
					if(current_edge_length < threshold_join_nodes &&
							  (time-final_time_nodes_created[current_node_global_ci])>protection_time &&
							  (time-final_time_nodes_created[clockwise_node_global_ci])>protection_time){
						
// 						printf("hello \n");
						
						// look for a cell that shares the edge - we need this for a T1 swap
						bool cell_with_same_edge_found = false;
						int cell_with_same_edge_mi = -1;
						
						unsigned counter_1 = 0;
						// loop over all cells containing current_node_global
						for(unsigned i=0;i<final_cells_per_node[current_node_global_ci];i++){
							
							unsigned temp_cell_1_mi = (unsigned)final_cell_store[current_node_global_ci + i*array_sizes];
							
							// if not current_cell, loop over all cells containing clockwise_node_global
							if(temp_cell_1_mi != current_cell_mi && !cell_with_same_edge_found){
								for(unsigned j=0;j<final_cells_per_node[clockwise_node_global_ci];j++){
									
									unsigned temp_cell_2_mi = (unsigned)final_cell_store[clockwise_node_global_ci + j*array_sizes];
									
									// if cell containing clockwise_node_global is also a cell containing current_node_global,
									// but not current_cell, then it is cell_with_same_edge;
									if(temp_cell_2_mi == temp_cell_1_mi){
										cell_with_same_edge_mi = temp_cell_1_mi;
										cell_with_same_edge_found = true;
										break;
									}
								}
							}
							
						}
						
						if(!cell_with_same_edge_found){
							
							join_logical = true;
							(*no_join_nodes_this_iteration)++;
							
// 							printf("doing a join \n");
							
							double new_node_position[2];
							for(unsigned dim=0;dim<2;dim++) {
								new_node_position[dim] = (current_node_position[dim] +
										  clockwise_node_position[dim])/2;
							}
							
							// find two unused nodes to put new nodes in
							unsigned new_node_ci = 0;
							while(final_cells_per_node[new_node_ci]>0){
								new_node_ci++;
							}
							unsigned new_node_mi = new_node_ci+1;
							
							assert(new_node_mi <= array_sizes);
							
							// store new node position
							for(unsigned dim=0;dim<2;dim++){
								final_node_positions[new_node_ci+dim*array_sizes] =
										  new_node_position[dim];
							}
							
							// set time new nodes created to current time
							final_time_nodes_created[new_node_ci] = time;
							
							// edit current cell nodes to contain new node instead of previous two nodes
							mxArray* mx_cell_nodes_edited = mxCreateDoubleMatrix(1, no_cell_nodes-1, mxREAL);
							double* cell_nodes_edited = mxGetPr(mx_cell_nodes_edited);
							
							int temp_counter = -1;
							for(unsigned temp_current_node_local=0;temp_current_node_local<no_cell_nodes;temp_current_node_local++){
								
								unsigned temp_current_node_global_mi = (unsigned)cell_nodes[temp_current_node_local];
								
								if(temp_current_node_global_mi==current_node_global_mi){
									temp_counter++;
									cell_nodes_edited[temp_counter]=new_node_mi;
								}
								else if(temp_current_node_global_mi==clockwise_node_global_mi){}
								else{
									temp_counter++;
									cell_nodes_edited[temp_counter]=temp_current_node_global_mi;
								}
							}
							mxSetCell(plhs[2], current_cell_ci, mx_cell_nodes_edited);
							
// 							printf("edited current cell \n");
							
							final_cells_per_node[new_node_ci] = 1;
							final_cell_store[new_node_ci] = current_cell_mi;
							
							if(final_cells_per_node[current_node_global_ci] > 1.1){
								
								for(unsigned i=0;i<final_cells_per_node[current_node_global_ci];i++){
									
									unsigned temp_cell_mi = (unsigned)final_cell_store[current_node_global_ci + i*array_sizes];
									
									if(temp_cell_mi!=current_cell_mi){
										
										unsigned temp_cell_ci = temp_cell_mi-1;
										
										mxArray* mx_cell_nodes_temp_cell = mxGetCell(plhs[2], temp_cell_ci);
										double* cell_nodes_temp_cell = mxGetPr(mx_cell_nodes_temp_cell);
										unsigned no_cell_nodes_temp_cell = mxGetN(mx_cell_nodes_temp_cell);
										
										mxArray* mx_cell_nodes_temp_cell_edited = mxCreateDoubleMatrix(1, no_cell_nodes_temp_cell, mxREAL);
										double* cell_nodes_temp_cell_edited = mxGetPr(mx_cell_nodes_temp_cell_edited);
										
										for(unsigned temp_current_node_local=0;temp_current_node_local<no_cell_nodes_temp_cell;temp_current_node_local++){
											
											unsigned temp_current_node_global_mi = (unsigned)cell_nodes_temp_cell[temp_current_node_local];
											
											if(temp_current_node_global_mi==current_node_global_mi){
												cell_nodes_temp_cell_edited[temp_current_node_local] = new_node_mi;
											}
											else{
												cell_nodes_temp_cell_edited[temp_current_node_local] = temp_current_node_global_mi;
											}
										}
										
										mxSetCell(plhs[2], temp_cell_ci, mx_cell_nodes_temp_cell_edited);
										
										final_cells_per_node[new_node_ci]++;
										unsigned matrix_index = (unsigned)(new_node_ci+array_sizes*(final_cells_per_node[new_node_ci]-1));
										final_cell_store[matrix_index] = temp_cell_mi;
									}
								}
// 								printf("edited cells with node 1 \n");
							}
							
							if(final_cells_per_node[clockwise_node_global_ci] > 1.1){
								
								for(unsigned i=0;i<final_cells_per_node[clockwise_node_global_ci];i++){
									
									unsigned temp_cell_mi = (unsigned)final_cell_store[clockwise_node_global_ci + i*array_sizes];
									
									if(temp_cell_mi!=current_cell_mi){
										
										unsigned temp_cell_ci = temp_cell_mi-1;
										
										// edit stretch cell 1 by adding in two new nodes in place of current node global
										mxArray* mx_cell_nodes_temp_cell = mxGetCell(plhs[2], temp_cell_ci);
										double* cell_nodes_temp_cell = mxGetPr(mx_cell_nodes_temp_cell);
										unsigned no_cell_nodes_temp_cell = mxGetN(mx_cell_nodes_temp_cell);
										
										mxArray* mx_cell_nodes_temp_cell_edited = mxCreateDoubleMatrix(1, no_cell_nodes_temp_cell, mxREAL);
										double* cell_nodes_temp_cell_edited = mxGetPr(mx_cell_nodes_temp_cell_edited);
										
										for(unsigned temp_current_node_local=0;temp_current_node_local<no_cell_nodes_temp_cell;temp_current_node_local++){
											
											unsigned temp_current_node_global_mi = (unsigned)cell_nodes_temp_cell[temp_current_node_local];
											
											if(temp_current_node_global_mi==clockwise_node_global_mi){
												cell_nodes_temp_cell_edited[temp_current_node_local] = new_node_mi;
											}
											else{
												cell_nodes_temp_cell_edited[temp_current_node_local] = temp_current_node_global_mi;
											}
										}
										
										mxSetCell(plhs[2], temp_cell_ci, mx_cell_nodes_temp_cell_edited);
										
										// add stretch cell 1 to cells_per_node and cell_store for new_node and new_node_2
										final_cells_per_node[new_node_ci]++;
										unsigned matrix_index = (unsigned)(new_node_ci+array_sizes*(final_cells_per_node[new_node_ci]-1));
										final_cell_store[matrix_index] = temp_cell_mi;
										
									}
								}
								printf("edited cells with node 2 \n");
								
							}
							
							// reset cells per node and cell store for current_node_global and clockwise_node_global.
							// they should no longer be part of any cells, unless something has gone horribly wrong.
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
						}
					}
				}
				if(join_logical){break;}
			}
		}
		if(join_logical){break;}
	}
}