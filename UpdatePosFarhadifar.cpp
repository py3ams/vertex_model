#include<cmath>
#include<cstdlib>
#include "mex.h"
#include "CommonFunctions.hpp"

#define pi 3.1415926535897932384626433832795

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if(nrhs != 11) {
        mexErrMsgTxt("11 input arguments required");
    }
    
    unsigned no_cells = mxGetM(prhs[0]);
    unsigned array_sizes = mxGetM(prhs[1]);
    
    double* initial_vertex_positions = mxGetPr(prhs[1]);
	double* boundary_element = mxGetPr(prhs[2]);
    double* cell_areas = mxGetPr(prhs[3]);
    double* cell_perimeters = mxGetPr(prhs[4]);
    double* cells_per_vertex = mxGetPr(prhs[5]);
    double delta_t = mxGetScalar(prhs[6]);
    double* target_areas = mxGetPr(prhs[7]);
    double* area_force_constants = mxGetPr(prhs[8]);
    double* perimeter_force_constants = mxGetPr(prhs[9]);
    double* tension_force_constants = mxGetPr(prhs[10]);
    
    plhs[0] = mxCreateDoubleMatrix(array_sizes, 2, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(array_sizes, 2, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    double* final_vertex_positions = mxGetPr(plhs[0]);
    double* vertex_movements = mxGetPr(plhs[1]);
    double* total_area_force = mxGetPr(plhs[2]);
    double* total_perimeter_force = mxGetPr(plhs[3]);
    double* total_tension_force = mxGetPr(plhs[4]);
    double* total_force = mxGetPr(plhs[5]);
    
    double area_force[array_sizes][2];
    double perimeter_force[array_sizes][2];
    double tension_force[array_sizes][2];
    double net_force[array_sizes][2];
    
    for(unsigned i=0;i<array_sizes;i++) {
        for(unsigned dim=0;dim<2;dim++) {
            area_force[i][dim] = 0;
            perimeter_force[i][dim] = 0;
            tension_force[i][dim] = 0;
            net_force[i][dim] = 0;
        }
    }
	
	// have moved this here temporarily. tension force constants should be related to 
	// junctions rather than cells.
	double area_force_constant = area_force_constants[0];
	double perimeter_force_constant = perimeter_force_constants[0];
	double tension_force_constant = tension_force_constants[0];
    
    // loop over all cells
    for(unsigned current_cell=0; current_cell<no_cells; current_cell++) {
        
        // get cell nodes
        mxArray* mx_cell_nodes = mxGetCell(prhs[0], current_cell);
        unsigned no_cell_nodes = mxGetN(mx_cell_nodes);
        
        if(no_cell_nodes>0){
            
            double* cell_nodes = mxGetPr(mx_cell_nodes);
            

            
// 		printf("start of main loop \n");
                       
            // get other stuff
            double current_cell_perimeter = cell_perimeters[current_cell];
            double current_cell_area = cell_areas[current_cell];
            
            // find height-to-area ratio and area^2 and average internal angle
            double current_cell_area_sq =  pow(current_cell_area, 2);
            
            double deviation_from_target_area = current_cell_area-target_areas[current_cell];
            
//		first loop over cell nodes. we will check for any concave nodes, and calculate the
//	   tension unit normals, during this loop.
            for(unsigned current_node_local=0;current_node_local<no_cell_nodes;current_node_local++) {
                
                bool concave_logical = false;
                
                // find the global value of the current cell node, and both its anti clockwise
                // and clockwise neighbours
                unsigned anti_clockwise_node_local = (current_node_local+no_cell_nodes-1)%no_cell_nodes;
                unsigned clockwise_node_local = (current_node_local+1)%no_cell_nodes;
                
                unsigned current_node_global_mi = (unsigned)cell_nodes[current_node_local];
                unsigned anti_clockwise_node_global_mi = (unsigned)cell_nodes[anti_clockwise_node_local];
                unsigned clockwise_node_global_mi = (unsigned)cell_nodes[clockwise_node_local];
                
                unsigned current_node_global_ci = current_node_global_mi-1;
                unsigned anti_clockwise_node_global_ci = anti_clockwise_node_global_mi-1;
                unsigned clockwise_node_global_ci = clockwise_node_global_mi-1;
                
                // extract the positions of the three nodes in question from the initial node positions matrix.
                double current_node_position[2];
                double anti_clockwise_node_position[2];
                double clockwise_node_position[2];
                
                for(unsigned dim=0;dim<2;dim++) {
                    
                    unsigned matrix_offset = dim*array_sizes;
                    
                    current_node_position[dim] = initial_vertex_positions[current_node_global_ci+matrix_offset];
                    anti_clockwise_node_position[dim] = initial_vertex_positions[anti_clockwise_node_global_ci+matrix_offset];
                    clockwise_node_position[dim] = initial_vertex_positions[clockwise_node_global_ci+matrix_offset];
                    
                }
                
                double tension_unit_normals_anti_clockwise[2];
                double tension_unit_normals_clockwise[2];
                
                // find unit vectors between the points
                findUnitVector(current_node_position, anti_clockwise_node_position,
                        tension_unit_normals_anti_clockwise);
                findUnitVector(current_node_position, clockwise_node_position,
                        tension_unit_normals_clockwise);
                
                // if the cross product is pointing in the opposite direction from the surface
                // unit normal, the current point is concave.
                if ((tension_unit_normals_anti_clockwise[0]*tension_unit_normals_clockwise[1] -
                        tension_unit_normals_clockwise[0]*tension_unit_normals_anti_clockwise[1]) < 0) {
                    
                    concave_logical = true;
                }
                
                double dot_product_tension_vectors = 0;
                double sum_tension_vectors[2];
                
                for(unsigned dim=0;dim<2;dim++) {
                    
                    sum_tension_vectors[dim] = tension_unit_normals_anti_clockwise[dim] +
                            tension_unit_normals_clockwise[dim];
                    
                    dot_product_tension_vectors += tension_unit_normals_anti_clockwise[dim]*
                            tension_unit_normals_clockwise[dim];
                    
                }
                
                // this can happen after mitosis at the new node - not really sure why
                if(dot_product_tension_vectors>1){
                    dot_product_tension_vectors = 1;
                }
                else if(dot_product_tension_vectors<-1){
                    dot_product_tension_vectors = -1;
                }

                double pressure_unit_normals[2];
                
// 			printf("%f \n", dot_product_tension_vectors);
// 			printf("%f \n", internal_angle);
                
                if(fabs(sum_tension_vectors[0])<0.000000001 && fabs(sum_tension_vectors[1])<0.000000001){
                    //this occurs just after mitosis - don't delete again!
//
                    pressure_unit_normals[0] = -tension_unit_normals_clockwise[1];
                    pressure_unit_normals[1] = tension_unit_normals_clockwise[0];
//
// 				printf("%f %f \n",pressure_unit_normals[0],pressure_unit_normals[1]);
//
                }
                
                else{
                    
                    double magnitude_resolved_tension = findNorm(sum_tension_vectors);
//				printf("%f \n",magnitude_resolved_tension);
                    
                    // if node is concave, pressure points in the same direction as the resolved tension. Otherwise,
                    // it points in the opposite direction.
                    if(concave_logical){
                        for(unsigned dim=0;dim<2;dim++) {
                            pressure_unit_normals[dim] =
                                    sum_tension_vectors[dim]/magnitude_resolved_tension;
                        }
                    }
                    else {
                        for(unsigned dim=0;dim<2;dim++) {
                            pressure_unit_normals[dim] =
                                    -sum_tension_vectors[dim]/magnitude_resolved_tension;
                        }
                    }
                }

                
                for(unsigned dim=0;dim<2;dim++) {
                                            
                    area_force[current_node_global_ci][dim] += -area_force_constant*
                            deviation_from_target_area*pressure_unit_normals[dim];
					
					tension_force[current_node_global_ci][dim] +=
							-tension_force_constant*tension_unit_normals_clockwise[dim];
					
                    perimeter_force[current_node_global_ci][dim] +=
                            -perimeter_force_constant*current_cell_perimeter*
                            (tension_unit_normals_clockwise[dim]+
                            tension_unit_normals_anti_clockwise[dim]);


                }
				
// 				if(current_cell==28){
// 					printf("%f %f \n", -tension_force_constant*tension_unit_normals_clockwise[0], -tension_force_constant*tension_unit_normals_clockwise[1]);
// 				}
            }
        }
    }
	
	/////////////////////////////////////////////////////////////////////////////////
	
	// we need to loop over the boundary nodes as some of the energy terms will not be 
	// found in the loop over cells. normally we can be sure that if a cell covers the ij
	// junction, the adjoining cell will contain ji. on the boundary there are no adjacent
	// cells, so we need to do these manually
	
	unsigned no_boundary_nodes = mxGetN(prhs[2]);
    
    for(unsigned current_node_local=0;current_node_local<no_boundary_nodes;current_node_local++){
        
// 		printf("hello 2 \n");
        
        bool concave_logical = false;
        
        unsigned anti_clockwise_node_local = (current_node_local+no_boundary_nodes-1)%no_boundary_nodes;
                
        unsigned current_node_global_mi = (unsigned)boundary_element[current_node_local];
        unsigned anti_clockwise_node_global_mi = (unsigned)boundary_element[anti_clockwise_node_local];
        
        unsigned current_node_global_ci = current_node_global_mi-1;
        unsigned anti_clockwise_node_global_ci = anti_clockwise_node_global_mi-1;
        
        double current_node_position[2];
        double anti_clockwise_node_position[2];
        
        for(unsigned dim=0;dim<2;dim++) {
            
            unsigned matrix_offset = dim*array_sizes;
            
            current_node_position[dim] = initial_vertex_positions[current_node_global_ci+matrix_offset];
            anti_clockwise_node_position[dim] = initial_vertex_positions[anti_clockwise_node_global_ci+matrix_offset];
            
        }
        
        double unit_vector_to_anti_clockwise_node[2];
        
        findUnitVector(current_node_position, anti_clockwise_node_position,
                unit_vector_to_anti_clockwise_node);
        
        for(unsigned dim=0;dim<2;dim++) {
            
			tension_force[current_node_global_ci][dim] +=
					-tension_force_constant*unit_vector_to_anti_clockwise_node[dim];
            
        }
    }
    

    /////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////// apply forces //////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
    
    *total_area_force = 0;
    *total_perimeter_force = 0;
    *total_tension_force = 0;
    *total_force = 0;
    
    for(unsigned current_node_global_ci=0;current_node_global_ci<array_sizes;current_node_global_ci++) {
        
        if(cells_per_vertex[current_node_global_ci]>0){
            
            double vertex_speed[2];
            
            for(unsigned dim=0;dim<2;dim++) {
                
                net_force[current_node_global_ci][dim] +=
                        area_force[current_node_global_ci][dim] +
                        perimeter_force[current_node_global_ci][dim] +
                        tension_force[current_node_global_ci][dim];
                
//                 if(current_node_global_ci==248){
//                     printf("%f ", area_force[current_node_global_ci][dim]);
//                 }
                
                unsigned matrix_index = current_node_global_ci+dim*array_sizes;
                
                vertex_movements[matrix_index] = net_force[current_node_global_ci][dim]*delta_t;
                
                final_vertex_positions[matrix_index] =
                        initial_vertex_positions[matrix_index] + vertex_movements[matrix_index];
                
            }
            
//             if(current_node_global_ci==248){
//                 printf("\n");
//             }
            
            double area_force_norm = findNorm(area_force[current_node_global_ci]);
            double perimeter_force_norm = findNorm(perimeter_force[current_node_global_ci]);
            double tension_force_norm = findNorm(tension_force[current_node_global_ci]);
            double net_force_norm = findNorm(net_force[current_node_global_ci]);
            
            *total_area_force += area_force_norm;
            *total_perimeter_force += perimeter_force_norm;
            *total_tension_force += tension_force_norm;
            *total_force += net_force_norm;
            
//             double vertex_speed_magnitude = findNorm(vertex_speed);
            
            
        }
    }
}
