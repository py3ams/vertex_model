#include<cmath>
#include<vector>
#include "mex.h"
#include "CommonFunctions.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double edge_threshold;

    if (nrhs < 3) {
        edge_threshold = 1e9;
    }
    else{
        edge_threshold = mxGetScalar(prhs[2]);
    }
    
    // find number of cells and size of node positions array
    unsigned length_cells = mxGetM(prhs[0]);
    unsigned array_sizes = mxGetM(prhs[1]);
    
    // get hold of the node positions
    double* node_positions = mxGetPr(prhs[1]);
    
    // create a matrix for each of the cell areas
    plhs[0] = mxCreateDoubleMatrix(length_cells, 1, mxREAL);
    double* cell_areas = mxGetPr(plhs[0]);
    
    // same for cell perimeters
    plhs[1] = mxCreateDoubleMatrix(length_cells, 1, mxREAL);
    double* cell_perimeters = mxGetPr(plhs[1]);
    
    // create cell matrix for each edge length
    plhs[2] = mxCreateCellMatrix(length_cells, 1);
    
    // create matrices for all the stats
    plhs[3] = mxCreateDoubleMatrix(1, 5, mxREAL);
    double* cell_area_stats = mxGetPr(plhs[3]);
    
    plhs[4] = mxCreateDoubleMatrix(1, 4, mxREAL);
    double* cell_perimeter_stats = mxGetPr(plhs[4]);
    
    plhs[5] = mxCreateDoubleMatrix(1, 4, mxREAL);
    double* shape_index_stats = mxGetPr(plhs[5]);
    
    plhs[6] = mxCreateDoubleMatrix(1, 3, mxREAL);
    double* edge_length_stats = mxGetPr(plhs[6]);
    
    plhs[7] = mxCreateDoubleMatrix(length_cells, 1, mxREAL);
    double* shape_indices = mxGetPr(plhs[7]);
    
    unsigned length_long_edges = 2*(array_sizes+length_cells);
    
    plhs[8] = mxCreateDoubleMatrix(length_long_edges, 2, mxREAL);
    double* long_edges = mxGetPr(plhs[8]);
        
    double max_edge_length = -100000;
    double min_edge_length = 100000;
    double mean_edge_length;
    
    double max_cell_perimeter = -100000;
    double min_cell_perimeter = 100000;
    double mean_cell_perimeter, std_cell_perimeter;
    
    double max_cell_area = -100000;
    double min_cell_area = 100000;
    double mean_cell_area, std_cell_area;
    
    double max_shape_index = -100000;
    double min_shape_index = 100000;
    double mean_shape_index, std_shape_index;
    
    double total_cell_perimeter = 0, total_cell_area = 0, total_shape_index = 0;
    double total_edge_length = 0;
    unsigned total_edges = 0, total_long_edges = 0;
    
    unsigned no_actual_cells = 0;
    
// loop over all the cells
    for(unsigned current_cell=0; current_cell<length_cells; current_cell++) {

        // find the current cell, and figure out the number of cell nodes and what
        // they are
        mxArray* mx_cell_nodes = mxGetCell(prhs[0], current_cell);
        unsigned no_cell_nodes = mxGetN(mx_cell_nodes);
        double* cell_nodes = mxGetPr(mx_cell_nodes);
        
        // create an mxArray to hold the edge lengths of this cell, and find the pointer to it
        mxArray* mx_current_cell_edge_lengths = mxCreateDoubleMatrix(no_cell_nodes, 1, mxREAL);
        double* current_cell_edge_lengths = mxGetPr(mx_current_cell_edge_lengths);
        
        double current_cell_perimeter = 0;
        
        // initialse array for the cell centre to 0
        double current_cell_centre[2] = {0};
        
        if(no_cell_nodes>0){
            
            no_actual_cells++;
            
            // first loop over the cell nodes. we will use this to find the edge lengths, perimeter,
            // and centre of the cell.
            for(unsigned current_node_local=0;current_node_local<no_cell_nodes;current_node_local++) {
                
                // find the current node and its clockwise counterpart in matlab indicies
                unsigned current_node_global_mi = (unsigned)cell_nodes[current_node_local];
                unsigned clockwise_node_global_mi = (unsigned)cell_nodes[(current_node_local+1)%no_cell_nodes];
                
                // convert to c indicies
                unsigned current_node_global_ci = current_node_global_mi-1;
                unsigned clockwise_node_global_ci = clockwise_node_global_mi-1;
                
                // extract the positions of the two nodes in question from the initial node positions matrix.
                double current_node_position[2], clockwise_node_position[2];
                
                for(unsigned dim=0;dim<2;dim++) {
                    
                    unsigned matrix_offset = dim*array_sizes;
                    
                    current_node_position[dim] = node_positions[current_node_global_ci + matrix_offset];
                    clockwise_node_position[dim] = node_positions[clockwise_node_global_ci + matrix_offset];
                }
                
                double current_edge_length =
                        findStraightLineDistanceBetweenTwoNodes(current_node_position, clockwise_node_position);
                
                if(current_edge_length>edge_threshold){
                    long_edges[total_long_edges] = current_node_global_mi;
                    long_edges[total_long_edges+length_long_edges] = clockwise_node_global_mi;
                    total_long_edges++;
                }                    
                
                // add the current edge length to the perimeter, and store the edge length in the
                // matrix of edge lengths for the current cell
                current_cell_perimeter += current_edge_length;
                current_cell_edge_lengths[current_node_local] = current_edge_length;
                
                max_edge_length = max(max_edge_length, current_edge_length);
                min_edge_length = min(min_edge_length, current_edge_length);
                
                total_edges++;
                total_edge_length += current_edge_length;
                
                // add the current node position to the centre variable (this is not the same as the centroid)
                for(unsigned dim=0;dim<2;dim++) {
                    current_cell_centre[dim] += current_node_position[dim]/(double)no_cell_nodes;
                }
            }
        }
        
// 		if(current_cell==length_cells-1){
// 			printf("%f %f \n",current_cell_centre[0],current_cell_centre[1]);
// 		}
        
        // add the matrix of edge lengths for the current cell into the matrix of cells storing
        // all the edge lengths
        mxSetCell(plhs[2], current_cell, mx_current_cell_edge_lengths);
        
        cell_perimeters[current_cell] = current_cell_perimeter;
        total_cell_perimeter += current_cell_perimeter;
        
        double current_cell_area = 0;
        
        if(no_cell_nodes>0){
            
            max_cell_perimeter = max(max_cell_perimeter, current_cell_perimeter);
            min_cell_perimeter = min(min_cell_perimeter, current_cell_perimeter);
                        
            // second loop over the cell nodes. we can now find the area of the cell (needed to find the
            // centre first)
            for(unsigned current_node_local=0;current_node_local<no_cell_nodes;current_node_local++) {
                
                // find the current node and its clockwise counterpart in matlab indicies
                unsigned current_node_global_mi = (unsigned)cell_nodes[current_node_local];
                unsigned clockwise_node_global_mi = (unsigned)cell_nodes[(current_node_local+1)%no_cell_nodes];
                
                // convert to c indicies
                unsigned current_node_global_ci = current_node_global_mi-1;
                unsigned clockwise_node_global_ci = clockwise_node_global_mi-1;
                
                // extract the positions of the two nodes in question from the initial node positions matrix.
                double current_node_position[3], clockwise_node_position[3];
                
                for(unsigned dim=0;dim<2;dim++) {
                    
                    unsigned matrix_offset = dim*array_sizes;
                    
                    current_node_position[dim] = node_positions[current_node_global_ci + matrix_offset];
                    clockwise_node_position[dim] = node_positions[clockwise_node_global_ci + matrix_offset];
                }
                
                // initialise some arrays we use to calculate the area
                double vector_current_node_to_centre[2];
                double vector_clockwise_node_to_centre[2];
                
                // find vectors from current node and clockwise node to centre of cell
                for(unsigned dim=0;dim<2;dim++) {
                    vector_current_node_to_centre[dim] = current_node_position[dim]-current_cell_centre[dim];
                    vector_clockwise_node_to_centre[dim] = clockwise_node_position[dim]-current_cell_centre[dim];
                }
                
                double a_sq = pow(vector_current_node_to_centre[0], 2) + pow(vector_current_node_to_centre[1], 2);
                double b_sq = pow(vector_clockwise_node_to_centre[0], 2) + pow(vector_clockwise_node_to_centre[1], 2);
                double c_sq = pow(current_cell_edge_lengths[current_node_local], 2);
                
                current_cell_area += 0.25*sqrt(4*a_sq*b_sq - pow(a_sq+b_sq-c_sq, 2));
                
// 			if(current_cell==1304){
// 				printf("%f %f %f \n", a_sq,b_sq,c_sq);
// 			}
                
            }
        }
        
        cell_areas[current_cell] = current_cell_area;
        total_cell_area += current_cell_area;
        
        // find the shape index (need to calculate both perimeter and area first)
        if(no_cell_nodes>0){
            max_cell_area = max(max_cell_area, current_cell_area);
            min_cell_area = min(min_cell_area, current_cell_area);
            double current_shape_index = pow(current_cell_perimeter, 2)/current_cell_area;
            shape_indices[current_cell] = current_shape_index;
            max_shape_index = max(max_shape_index, current_shape_index);
            min_shape_index = min(min_shape_index, current_shape_index);
            total_shape_index += current_shape_index;
        }
        else{
            shape_indices[current_cell] = 0;
        }
        
    }
    
    mean_cell_area = total_cell_area/(double)no_actual_cells;
    mean_cell_perimeter = total_cell_perimeter/(double)no_actual_cells;
    mean_edge_length = total_edge_length/(double)total_edges;
    mean_shape_index = total_shape_index/(double)no_actual_cells;
    
    if(nlhs>3){
        
        double std_cell_areas_sum = 0;
        double std_cell_perimeters_sum = 0;
        double std_shape_indices_sum = 0;
        
        for(unsigned current_cell_ci=0;current_cell_ci<length_cells;current_cell_ci++){
            
            if(cell_areas[current_cell_ci]>0){
                
                std_cell_areas_sum += pow(cell_areas[current_cell_ci]-mean_cell_area, 2);
                std_cell_perimeters_sum += pow(cell_perimeters[current_cell_ci]-mean_cell_perimeter, 2);
                std_shape_indices_sum += pow(shape_indices[current_cell_ci]-mean_shape_index, 2);
                
            }
            
        }
        
        std_cell_area = sqrt(std_cell_areas_sum/(double)no_actual_cells);
        std_cell_perimeter = sqrt(std_cell_perimeters_sum/(double)no_actual_cells);
        std_shape_index = sqrt(std_shape_indices_sum/(double)no_actual_cells);
        
        cell_perimeter_stats[0] = mean_cell_perimeter;
        cell_perimeter_stats[1] = max_cell_perimeter;
        cell_perimeter_stats[2] = min_cell_perimeter;
        cell_perimeter_stats[3] = std_cell_perimeter;
        
        cell_area_stats[0] = mean_cell_area;
        cell_area_stats[1] = max_cell_area;
        cell_area_stats[2] = min_cell_area;
        cell_area_stats[3] = std_cell_area;
        cell_area_stats[4] = total_cell_area;
        
        shape_index_stats[0] = mean_shape_index;
        shape_index_stats[1] = max_shape_index;
        shape_index_stats[2] = min_shape_index;
        shape_index_stats[3] = std_shape_index;
        
        edge_length_stats[0] = mean_edge_length;
        edge_length_stats[1] = max_edge_length;
        edge_length_stats[2] = min_edge_length;
        
    }
}
