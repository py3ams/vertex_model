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
    
    //	Get current_radius
    const mxArray* mx_current_radius = mexGetVariablePtr("global", "current_radius");
    if (mx_current_radius == NULL) {
        mexErrMsgTxt("current_radius not defined");
    }
    double current_radius = mxGetScalar(mx_current_radius);
    
    // Get elongation //
    const mxArray* mx_elongation = mexGetVariablePtr("global", "elongation");
    if (mx_elongation == NULL) {
        mexErrMsgTxt("elongation not defined");
    }
    double elongation = mxGetScalar(mx_elongation);
    
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
    
    //	Get split_probability
    const mxArray* mx_split_probability = mexGetVariablePtr("caller", "split_probability");
    if (mx_split_probability == NULL) {
        mexErrMsgTxt("split_probability not defined");
    }
    double split_probability = mxGetScalar(mx_split_probability);
    
    //	Get threshold_join_nodes
    const mxArray* mx_threshold_join_nodes = mexGetVariablePtr("caller", "threshold_join_nodes");
    if (mx_threshold_join_nodes == NULL) {
        mexErrMsgTxt("threshold_join_nodes not defined");
    }
    double threshold_join_nodes = mxGetScalar(mx_threshold_join_nodes);
    
    //	Get threshold_T1_swaps
    const mxArray* mx_threshold_T1_swaps = mexGetVariablePtr("caller", "threshold_T1_swaps");
    if (mx_threshold_T1_swaps == NULL) {
        mexErrMsgTxt("threshold_T1_swaps not defined");
    }
    double threshold_T1_swaps = mxGetScalar(mx_threshold_T1_swaps);
    
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
    plhs[4] = mxCreateDoubleMatrix(array_sizes, 3, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(array_sizes, 1, mxREAL);
    
    double* final_cell_store = mxGetPr(plhs[0]);
    double* final_cells_per_node = mxGetPr(plhs[2]);
    double* no_split_nodes_this_iteration = mxGetPr(plhs[3]);
    double* final_node_positions = mxGetPr(plhs[4]);
    double* final_time_nodes_created = mxGetPr(plhs[5]);
    
    *no_split_nodes_this_iteration = 0;
    
    // set output matrices equal to input matrices
    for(unsigned current_node=0;current_node<array_sizes;current_node++){
        final_cells_per_node[current_node] = initial_cells_per_node[current_node];
        final_time_nodes_created[current_node] = initial_time_nodes_created[current_node];
        for(unsigned dim=0;dim<3;dim++){
            final_node_positions[current_node+array_sizes*dim] = initial_node_positions[current_node+array_sizes*dim];
        }
        for(unsigned i=0;i<width_cell_store;i++){
            final_cell_store[current_node+i*array_sizes] = initial_cell_store[current_node+i*array_sizes];
        }
    }
    
    // set output cells equal to input cells
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
    
    // main_loop
    for(unsigned current_node_global_ci=0;current_node_global_ci<array_sizes;current_node_global_ci++){
        
        unsigned current_node_global_mi = current_node_global_ci+1;
        
        double rand_number = rand()/((double)RAND_MAX);
        unsigned no_cells_with_current_node = (unsigned)final_cells_per_node[current_node_global_ci];
        
        if((no_cells_with_current_node>=5) && (time-final_time_nodes_created[current_node_global_ci]>protection_time) &&
                (rand_number < split_probability)){
            
            (*no_split_nodes_this_iteration)++;
            
            double current_node_position[3];
            
            for(unsigned dim=0;dim<3;dim++) {
                current_node_position[dim] = final_node_positions[current_node_global_ci + dim*array_sizes];
            }
            
            double surface_unit_normal[3];
            findSurfaceUnitNormal(current_node_position, current_radius, elongation, surface_unit_normal);
            
            unsigned stretch_cell_1_local = rand()%no_cells_with_current_node;
            unsigned stretch_cell_1_global_mi = (unsigned)final_cell_store[current_node_global_ci+array_sizes*stretch_cell_1_local];
            unsigned stretch_cell_1_global_ci = stretch_cell_1_global_mi-1;
            
            mxArray* mx_cell_nodes_stretch_cell_1 = mxGetCell(plhs[1], stretch_cell_1_global_ci);
            double* cell_nodes_stretch_cell_1 = mxGetPr(mx_cell_nodes_stretch_cell_1);
            unsigned no_cell_nodes_stretch_cell_1 = mxGetN(mx_cell_nodes_stretch_cell_1);
            
            unsigned possible_stretch_cells_2[no_cells_with_current_node-1];
            unsigned no_possible_stretch_cells_2 = 0;
            
            for(unsigned i=0; i<no_cells_with_current_node; i++){
                
                unsigned temp_cell_mi = (unsigned)final_cell_store[current_node_global_ci+i*array_sizes];
                unsigned temp_cell_ci = temp_cell_mi-1;
                
                mxArray* mx_temp_cell_nodes = mxGetCell(plhs[1], temp_cell_ci);
                double* temp_cell_nodes = mxGetPr(mx_temp_cell_nodes);
                unsigned no_temp_cell_nodes = mxGetN(mx_temp_cell_nodes);
                
                unsigned no_common_nodes = 0;
                
                for(unsigned j=0;j<no_temp_cell_nodes;j++){
                    for(unsigned k=0;k<no_cell_nodes_stretch_cell_1;k++){
                        if(temp_cell_nodes[j]==cell_nodes_stretch_cell_1[k]){   
                            no_common_nodes++;
                        }
                    }
                }
                
                if(no_common_nodes==1){
                    no_possible_stretch_cells_2++;
                    possible_stretch_cells_2[no_possible_stretch_cells_2-1] = temp_cell_mi;
                }
            }
            
            unsigned stretch_cell_2_local = rand()%no_possible_stretch_cells_2;
            unsigned stretch_cell_2_global_mi = possible_stretch_cells_2[stretch_cell_2_local];
            unsigned stretch_cell_2_global_ci = stretch_cell_2_global_mi-1;
            
            mxArray* mx_cell_nodes_stretch_cell_2 = mxGetCell(plhs[1], stretch_cell_2_global_ci);
            double* cell_nodes_stretch_cell_2 = mxGetPr(mx_cell_nodes_stretch_cell_2);
            unsigned no_cell_nodes_stretch_cell_2 = mxGetN(mx_cell_nodes_stretch_cell_2);
            
            double centre_stretch_cell_1[3] = {0};
                        
            for(unsigned i=0;i<no_cell_nodes_stretch_cell_1;i++){
                
                for(unsigned dim=0;dim<3;dim++){
                    unsigned matrix_index = (unsigned)cell_nodes_stretch_cell_1[i]-1+dim*array_sizes;
                    centre_stretch_cell_1[dim] += 
                            final_node_positions[matrix_index]/no_cell_nodes_stretch_cell_1;
                }
            }
            
            double centre_stretch_cell_2[3] = {0};
            
            for(unsigned i=0;i<no_cell_nodes_stretch_cell_2;i++){
                
                for(unsigned dim=0;dim<3;dim++){
                    unsigned matrix_index = (unsigned)cell_nodes_stretch_cell_2[i]-1+dim*array_sizes;
                    centre_stretch_cell_2[dim] +=
                            final_node_positions[matrix_index]/no_cell_nodes_stretch_cell_2;
                }
            }
            
            double mapped_centre_stretch_cell_1[3], mapped_centre_stretch_cell_2[3];
            
            map_node_position(centre_stretch_cell_1, current_radius, elongation,
                    mapped_centre_stretch_cell_1);
            map_node_position(centre_stretch_cell_2, current_radius, elongation,
                    mapped_centre_stretch_cell_2);
            
            double distance_stretch_cell_1_centre_to_plane = 0;
            double distance_stretch_cell_2_centre_to_plane = 0;
            
            for(unsigned dim=0;dim<3;dim++) {
                distance_stretch_cell_1_centre_to_plane += (current_node_position[dim] -
                        mapped_centre_stretch_cell_1[dim])*surface_unit_normal[dim];
                
                distance_stretch_cell_2_centre_to_plane += (current_node_position[dim] -
                        mapped_centre_stretch_cell_2[dim])*surface_unit_normal[dim];
            }
            
            double stretch_cell_1_centre_projection_to_plane[3], stretch_cell_2_centre_projection_to_plane[3];
            for(unsigned dim=0;dim<3;dim++) {
                stretch_cell_1_centre_projection_to_plane[dim] = mapped_centre_stretch_cell_1[dim] +
                        distance_stretch_cell_1_centre_to_plane*surface_unit_normal[dim];
                stretch_cell_2_centre_projection_to_plane[dim] = mapped_centre_stretch_cell_2[dim] +
                        distance_stretch_cell_2_centre_to_plane*surface_unit_normal[dim];
            }
            
            double unit_vector_on_plane[3];
                            
            findUnitVector(stretch_cell_1_centre_projection_to_plane, stretch_cell_2_centre_projection_to_plane,
                    unit_vector_on_plane);
            
            double new_node_1_direction[3];
            findCrossProduct(unit_vector_on_plane,
                    surface_unit_normal, new_node_1_direction);
            
            double norm_direction_vector = findNorm(new_node_1_direction);
                            
            double new_node_1_position_on_plane[3], new_node_2_position_on_plane[3];
            for(unsigned dim=0;dim<3;dim++){
                new_node_1_position_on_plane[dim] = current_node_position[dim] +
                        0.5*threshold_T1_swaps*new_node_1_direction[dim]/norm_direction_vector;
                new_node_2_position_on_plane[dim] = current_node_position[dim] -
                        0.5*threshold_T1_swaps*new_node_1_direction[dim]/norm_direction_vector;
            }
            
            double new_node_1_position_mapped[3], new_node_2_position_mapped[3];
            
            // map new node positions back to ellipsoid surface
            map_node_position(new_node_1_position_on_plane, current_radius, elongation,
                    new_node_1_position_mapped);
            map_node_position(new_node_2_position_on_plane, current_radius, elongation,
                    new_node_2_position_mapped);
            
            // find two unused nodes to put new nodes in
            unsigned new_node_1_ci = 0;
            while(final_cells_per_node[new_node_1_ci]>0){
                new_node_1_ci++;
            }
            unsigned new_node_1_mi = new_node_1_ci+1;
            
            unsigned new_node_2_ci = new_node_1_ci+1;
            while(final_cells_per_node[new_node_2_ci]>0){
                new_node_2_ci++;
            }
            unsigned new_node_2_mi = new_node_2_ci+1;
            
            // store new node positions
            for(unsigned dim=0;dim<3;dim++){
                final_node_positions[new_node_1_ci+dim*array_sizes] =
                        new_node_1_position_mapped[dim];
                final_node_positions[new_node_2_ci+dim*array_sizes] =
                        new_node_2_position_mapped[dim];
            }

                
            mxArray* mx_cell_nodes_stretch_cell_1_edited = mxCreateDoubleMatrix(1, no_cell_nodes_stretch_cell_1+1, mxREAL);
            double* cell_nodes_stretch_cell_1_edited = mxGetPr(mx_cell_nodes_stretch_cell_1_edited);
            
            unsigned index_new_node_1_in_stretch_cell_1_ci, index_new_node_2_in_stretch_cell_1_ci;
            int temp_counter = -1;
            for(unsigned temp_current_node_local=0;temp_current_node_local<no_cell_nodes_stretch_cell_1;temp_current_node_local++){
                
                temp_counter++;
                unsigned temp_current_node_global_mi = (unsigned)cell_nodes_stretch_cell_1[temp_current_node_local];
                
                if(temp_current_node_global_mi==current_node_global_mi){
                    cell_nodes_stretch_cell_1_edited[temp_counter] = new_node_2_mi;
                    index_new_node_2_in_stretch_cell_1_ci = temp_counter;
                    temp_counter++;
                    cell_nodes_stretch_cell_1_edited[temp_counter] = new_node_1_mi;
                    index_new_node_1_in_stretch_cell_1_ci = temp_counter;
                }
                else{
                    cell_nodes_stretch_cell_1_edited[temp_counter] = temp_current_node_global_mi;
                }
//                                 printf("%u \n",temp_counter);
            }
            
            mxSetCell(plhs[1], stretch_cell_1_global_ci, mx_cell_nodes_stretch_cell_1_edited);
            no_cell_nodes_stretch_cell_1++;
            
            final_cells_per_node[new_node_1_ci]++;
            unsigned matrix_index = (unsigned)(new_node_1_ci+array_sizes*(final_cells_per_node[new_node_1_ci]-1));
            final_cell_store[matrix_index] = stretch_cell_1_global_mi;
            final_cells_per_node[new_node_2_ci]++;
            matrix_index = (unsigned)(new_node_2_ci+array_sizes*(final_cells_per_node[new_node_2_ci]-1));
            final_cell_store[matrix_index] = stretch_cell_1_global_mi;
            
            // edit stretch cell 2 by adding in two new nodes in place of clockwise node global
            mxArray* mx_cell_nodes_stretch_cell_2_edited = mxCreateDoubleMatrix(1, no_cell_nodes_stretch_cell_2+1, mxREAL);
            double* cell_nodes_stretch_cell_2_edited = mxGetPr(mx_cell_nodes_stretch_cell_2_edited);
            
            unsigned index_new_node_1_in_stretch_cell_2_ci, index_new_node_2_in_stretch_cell_2_ci;
            temp_counter = -1;
            for(unsigned temp_current_node_local=0;temp_current_node_local<no_cell_nodes_stretch_cell_2;temp_current_node_local++){
                
                temp_counter++;
                unsigned temp_current_node_global_mi = (unsigned)cell_nodes_stretch_cell_2[temp_current_node_local];
                
//                                 printf("%u \n",temp_current_node_global_mi);
                
                if(temp_current_node_global_mi==current_node_global_mi){
                    cell_nodes_stretch_cell_2_edited[temp_counter] = new_node_1_mi;
                    index_new_node_1_in_stretch_cell_2_ci = temp_counter;
                    temp_counter++;
                    cell_nodes_stretch_cell_2_edited[temp_counter] = new_node_2_mi;
                    index_new_node_2_in_stretch_cell_2_ci = temp_counter;
                }
                else{
                    cell_nodes_stretch_cell_2_edited[temp_counter] = temp_current_node_global_mi;
                }
            }
            mxSetCell(plhs[1], stretch_cell_2_global_ci, mx_cell_nodes_stretch_cell_2_edited);
            no_cell_nodes_stretch_cell_2++;
            
//                             printf("Edited stretch cells \n");
            
            // add stretch cell 2 to cells_per_node and cell_store for new_node_1 and new_node_2
            final_cells_per_node[new_node_1_ci]++;
            matrix_index = (unsigned)(new_node_1_ci+array_sizes*(final_cells_per_node[new_node_1_ci]-1));
            final_cell_store[matrix_index] = stretch_cell_2_global_mi;
            final_cells_per_node[new_node_2_ci]++;
            matrix_index = (unsigned)(new_node_2_ci+array_sizes*(final_cells_per_node[new_node_2_ci]-1));
            final_cell_store[matrix_index] = stretch_cell_2_global_mi;
            
            unsigned temp_clockwise_node_local, temp_clockwise_node_global_mi;
            unsigned temp_anti_clockwise_node_local, temp_anti_clockwise_node_global_mi;
             
            // replace current_node_global with new node 1 in cells anti clockwise to stretch cell 1
            bool temp_flag = true;
            
            unsigned last_cell_edited_mi = stretch_cell_1_global_mi;
            
            if(index_new_node_1_in_stretch_cell_1_ci==no_cell_nodes_stretch_cell_1-1) {
                temp_clockwise_node_local = 0;
            }
            else {
                temp_clockwise_node_local = index_new_node_1_in_stretch_cell_1_ci+1;
            }
            
            temp_clockwise_node_global_mi = (unsigned)cell_nodes_stretch_cell_1_edited[temp_clockwise_node_local];
            
            while(temp_flag){
                
                temp_flag = false;
                
                for(unsigned temp_cell_local=0;temp_cell_local<no_cells_with_current_node;temp_cell_local++){
                    
                    unsigned matrix_index = (unsigned)current_node_global_ci+temp_cell_local*array_sizes;
                    
                    unsigned temp_cell_global_mi = (unsigned)final_cell_store[matrix_index];
                    unsigned temp_cell_global_ci = temp_cell_global_mi-1;
                    
                    if(temp_cell_global_mi!=last_cell_edited_mi && temp_cell_global_mi!=stretch_cell_2_global_mi){
                        
                        mxArray* mx_temp_cell_nodes = mxGetCell(plhs[1], temp_cell_global_ci);
                        double* temp_cell_nodes = mxGetPr(mx_temp_cell_nodes);
                        unsigned no_temp_cell_nodes = mxGetN(mx_temp_cell_nodes);
                        
                        for(unsigned temp_cell_node_local=0;temp_cell_node_local<no_temp_cell_nodes;temp_cell_node_local++){
                            if(temp_cell_nodes[temp_cell_node_local]==temp_clockwise_node_global_mi){
                                temp_flag = true;
                                break;
                            }
                        }
                    
                        if(temp_flag){
                            
                            mxArray* mx_temp_cell_nodes_edited = mxCreateDoubleMatrix(1, no_temp_cell_nodes, mxREAL);
                            double* temp_cell_nodes_edited = mxGetPr(mx_temp_cell_nodes_edited);
                            
                            unsigned index_new_node_1_in_temp_cell_ci;
                            for(unsigned temp_cell_node_local=0;temp_cell_node_local<no_temp_cell_nodes;temp_cell_node_local++){
                                if(temp_cell_nodes[temp_cell_node_local]==current_node_global_mi){
                                    temp_cell_nodes_edited[temp_cell_node_local] = new_node_1_mi;
                                    index_new_node_1_in_temp_cell_ci = temp_cell_node_local;
                                }
                                else{
                                    temp_cell_nodes_edited[temp_cell_node_local] = temp_cell_nodes[temp_cell_node_local];
                                }
                            }
                            mxSetCell(plhs[1], temp_cell_global_ci, mx_temp_cell_nodes_edited);
                            
                            final_cells_per_node[new_node_1_ci]++;
                            matrix_index = (unsigned)(new_node_1_ci+array_sizes*(final_cells_per_node[new_node_1_ci]-1));
                            final_cell_store[matrix_index] = temp_cell_global_mi;
                            
                            if(index_new_node_1_in_temp_cell_ci==(no_temp_cell_nodes-1)) {
                                temp_clockwise_node_local = 0;
                            }
                            else {
                                temp_clockwise_node_local = index_new_node_1_in_temp_cell_ci+1;
                            }
                            
                            temp_clockwise_node_global_mi = (unsigned)temp_cell_nodes[temp_clockwise_node_local];
                            last_cell_edited_mi = temp_cell_global_mi;
                            break;
                        }
                    }
                }
            }

            
            // replace current_node_global with new node 2 in cells clockwise to stretch cell 1
            temp_flag = true;
            
            last_cell_edited_mi = stretch_cell_1_global_mi;
            
            if(index_new_node_2_in_stretch_cell_1_ci == 0){
                temp_anti_clockwise_node_local = no_cell_nodes_stretch_cell_1-1;
            }
            else{
                temp_anti_clockwise_node_local = index_new_node_2_in_stretch_cell_1_ci-1;
            }
            
            temp_anti_clockwise_node_global_mi = (unsigned)cell_nodes_stretch_cell_1_edited[temp_anti_clockwise_node_local];
            
            while(temp_flag){
                
                temp_flag = false;
                
                for(unsigned temp_cell_local=0;temp_cell_local<no_cells_with_current_node;temp_cell_local++){
                    
                    unsigned matrix_index = (unsigned)current_node_global_ci+temp_cell_local*array_sizes;
                    
                    unsigned temp_cell_global_mi = (unsigned)final_cell_store[matrix_index];
                    unsigned temp_cell_global_ci = temp_cell_global_mi-1;
                    
                    if(temp_cell_global_mi!=last_cell_edited_mi && temp_cell_global_mi!=stretch_cell_2_global_mi){
                                           
                        mxArray* mx_temp_cell_nodes = mxGetCell(plhs[1], temp_cell_global_ci);
                        double* temp_cell_nodes = mxGetPr(mx_temp_cell_nodes);
                        unsigned no_temp_cell_nodes = mxGetN(mx_temp_cell_nodes);
                        
                        for(unsigned temp_cell_node_local=0;temp_cell_node_local<no_temp_cell_nodes;temp_cell_node_local++){
                            if(temp_cell_nodes[temp_cell_node_local]==temp_anti_clockwise_node_global_mi){
                                temp_flag = true;
                                break;
                            }
                        }

                        if(temp_flag){
                            
                            mxArray* mx_temp_cell_nodes_edited = mxCreateDoubleMatrix(1, no_temp_cell_nodes, mxREAL);
                            double* temp_cell_nodes_edited = mxGetPr(mx_temp_cell_nodes_edited);
                            
                            unsigned index_new_node_2_in_temp_cell_ci;
                            for(unsigned temp_cell_node_local=0;temp_cell_node_local<no_temp_cell_nodes;temp_cell_node_local++){
                                if(temp_cell_nodes[temp_cell_node_local]==current_node_global_mi){
                                    temp_cell_nodes_edited[temp_cell_node_local] = new_node_2_mi;
                                    index_new_node_2_in_temp_cell_ci = temp_cell_node_local;
                                }
                                else{
                                    temp_cell_nodes_edited[temp_cell_node_local] = temp_cell_nodes[temp_cell_node_local];
                                }
                            }
                            mxSetCell(plhs[1], temp_cell_global_ci, mx_temp_cell_nodes_edited);
                                                        
                            final_cells_per_node[new_node_2_ci]++;
                            matrix_index = (unsigned)(new_node_2_ci+array_sizes*(final_cells_per_node[new_node_2_ci]-1));
                            final_cell_store[matrix_index] = temp_cell_global_mi;
                            
                            if(index_new_node_2_in_temp_cell_ci==0) {
                                temp_anti_clockwise_node_local = no_temp_cell_nodes-1;
                            }
                            else {
                                temp_anti_clockwise_node_local = index_new_node_2_in_temp_cell_ci-1;
                            }
                            
                            temp_anti_clockwise_node_global_mi = (unsigned)temp_cell_nodes[temp_anti_clockwise_node_local];
                            last_cell_edited_mi = temp_cell_global_mi;
                            break;
                        }
                    }
                }
            }
            
            final_cells_per_node[current_node_global_ci] = 0;
            
            for(unsigned i=0;i<width_cell_store;i++){
                final_cell_store[current_node_global_ci+i*array_sizes] = 0;
            }
            
            for(unsigned dim=0;dim<3;dim++){
                final_node_positions[current_node_global_ci+dim*array_sizes] = 0;
            }
        }
    }
}

