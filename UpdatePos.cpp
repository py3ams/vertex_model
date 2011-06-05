#include<cmath>
#include<cstdlib>
#include "mex.h"
#include "CommonFunctions.hpp"

#define pi 3.1415926535897932384626433832795

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if(nrhs != 20) {
        mexErrMsgTxt("20 input arguments required");
    }
    
    unsigned no_cells = mxGetM(prhs[0]);
    unsigned array_sizes = mxGetM(prhs[1]);
    
    double* initial_node_positions = mxGetPr(prhs[1]);
    double* boundary_element = mxGetPr(prhs[2]);
    double* cell_areas = mxGetPr(prhs[3]);
    double* cell_perimeters = mxGetPr(prhs[4]);
    double* cell_volumes = mxGetPr(prhs[5]);
    double* cells_per_node = mxGetPr(prhs[6]);
    double delta_t = mxGetScalar(prhs[7]);
    //edge_lengths are prhs[8]
    double mean_edge_length = mxGetScalar(prhs[9]);
    double* target_areas = mxGetPr(prhs[10]);
    double tension_anisotropy_factor = mxGetScalar(prhs[11]);
    double viscosity = mxGetScalar(prhs[12]);
    double* area_force_constants = mxGetPr(prhs[13]);
    double boundary_deformation_force_constant = mxGetScalar(prhs[14]);
    double boundary_edge_force_constant = mxGetScalar(prhs[15]);
    double* deformation_force_constants = mxGetPr(prhs[16]);
    double* elongation_force_constants = mxGetPr(prhs[17]);
    double* perimeter_force_constants = mxGetPr(prhs[18]);
    double* tension_force_constants = mxGetPr(prhs[19]);
    
    unsigned no_internal_angles = 0;
    
    for(unsigned current_cell=0; current_cell<no_cells; current_cell++) {
        mxArray* mx_cell_nodes = mxGetCell(prhs[0], current_cell);
        no_internal_angles += mxGetN(mx_cell_nodes);
    }
    
    plhs[0] = mxCreateDoubleMatrix(array_sizes, 2, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(no_internal_angles, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(array_sizes, 2, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[6] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[7] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[8] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[9] = mxCreateDoubleMatrix(1, 1, mxREAL);
    plhs[10] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    double* final_node_positions = mxGetPr(plhs[0]);
    double* angle_deviations = mxGetPr(plhs[1]);
    double* vertex_movements = mxGetPr(plhs[2]);
    double* total_area_force = mxGetPr(plhs[3]);
    double* total_boundary_deformation_force = mxGetPr(plhs[4]);
    double* total_boundary_edge_force = mxGetPr(plhs[5]);
    double* total_deformation_force = mxGetPr(plhs[6]);
    double* total_elongation_force = mxGetPr(plhs[7]);
    double* total_perimeter_force = mxGetPr(plhs[8]);
    double* total_tension_force = mxGetPr(plhs[9]);
    double* total_force = mxGetPr(plhs[10]);
    
    double area_force[array_sizes][2];
    double boundary_deformation_force[array_sizes][2];
    double boundary_edge_force[array_sizes][2];
    double deformation_force[array_sizes][2];
    double elongation_force[array_sizes][2];
    double perimeter_force[array_sizes][2];
    double tension_force[array_sizes][2];
    double net_force[array_sizes][2];
    
    for(unsigned i=0;i<array_sizes;i++) {
        for(unsigned dim=0;dim<2;dim++) {
            area_force[i][dim] = 0;
            boundary_deformation_force[i][dim] = 0;
            boundary_edge_force[i][dim] = 0;
            deformation_force[i][dim] = 0;
            elongation_force[i][dim] = 0;
            perimeter_force[i][dim] = 0;
            tension_force[i][dim] = 0;
            net_force[i][dim] = 0;
        }
    }
    
    unsigned angle_counter = 0;
    
    // loop over all cells
    for(unsigned current_cell=0; current_cell<no_cells; current_cell++) {
        
        // get cell nodes
        mxArray* mx_cell_nodes = mxGetCell(prhs[0], current_cell);
        unsigned no_cell_nodes = mxGetN(mx_cell_nodes);
        
        if(no_cell_nodes>0){
            
            double* cell_nodes = mxGetPr(mx_cell_nodes);
            
            double area_force_constant = area_force_constants[current_cell];
            double deformation_force_constant = deformation_force_constants[current_cell];
            double elongation_force_constant = elongation_force_constants[current_cell];
            double perimeter_force_constant = perimeter_force_constants[current_cell];
            double tension_force_constant = tension_force_constants[current_cell];
            
// 		printf("start of main loop \n");
            
            // get edge lengths
            mxArray* mx_edge_lengths = mxGetCell(prhs[8], current_cell);
            double* edge_lengths = mxGetPr(mx_edge_lengths);
            
            // get other stuff
            double current_cell_perimeter = cell_perimeters[current_cell];
            double current_cell_area = cell_areas[current_cell];
            double current_cell_height =
                    cell_volumes[current_cell]/cell_areas[current_cell];
            
            // find height-to-area ratio and area^2 and average internal angle
            double current_cell_height_area_ratio = current_cell_height/current_cell_area;
            double current_cell_area_sq =  pow(current_cell_area, 2);
            double target_internal_angle = pi*(no_cell_nodes-2)/no_cell_nodes;
            
            double deviation_from_target_area = target_areas[current_cell]-current_cell_area;
            
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
                    
                    current_node_position[dim] = initial_node_positions[current_node_global_ci+matrix_offset];
                    anti_clockwise_node_position[dim] = initial_node_positions[anti_clockwise_node_global_ci+matrix_offset];
                    clockwise_node_position[dim] = initial_node_positions[clockwise_node_global_ci+matrix_offset];
                    
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
                
                double internal_angle;
                
                if(concave_logical){
                    internal_angle = 2*pi-acos(dot_product_tension_vectors);
                }
                else {
                    internal_angle = acos(dot_product_tension_vectors);
                }
                
                double deviation_from_target_angle = internal_angle-target_internal_angle;
                
                angle_deviations[angle_counter] = fabs(deviation_from_target_angle);
                angle_counter++;
                
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
                
                double clockwise_edge_angle_from_horizontal =
                        atan(fabs((current_node_position[1]-clockwise_node_position[1])/
                        (current_node_position[0]-clockwise_node_position[0])));
                
                double anti_clockwise_edge_angle_from_horizontal =
                        atan(fabs((current_node_position[1]-anti_clockwise_node_position[1])/
                        (current_node_position[0]-anti_clockwise_node_position[0])));
                
                double clockwise_edge_tension_factor =
                        tension_force_constant*(1+clockwise_edge_angle_from_horizontal*2/pi*tension_anisotropy_factor);
                
                double anti_clockwise_edge_tension_factor =
                        tension_force_constant*(1+anti_clockwise_edge_angle_from_horizontal*2/pi*tension_anisotropy_factor);
                
//             printf("%f %f \n",1+clockwise_edge_angle_from_horizontal*2/pi*tension_anisotropy_factor,clockwise_edge_tension_factor);
                
                for(unsigned dim=0;dim<2;dim++) {
                    
                    if(fabs(deviation_from_target_area)>1e-9){
                    
                        area_force[current_node_global_ci][dim] += area_force_constant*fabs(pow(
                                deviation_from_target_area, 3))/deviation_from_target_area*
                                pressure_unit_normals[dim];
                        
                    }
                    
                    if(fabs(deviation_from_target_angle)>1e-9){

                        deformation_force[current_node_global_ci][dim] += deformation_force_constant*
                                fabs(pow(deviation_from_target_angle, 3))/deviation_from_target_angle*
                                pressure_unit_normals[dim];
                        
                    }
                    
                    elongation_force[current_node_global_ci][dim] += elongation_force_constant*
                            pow(current_cell_height_area_ratio, 1)*pressure_unit_normals[dim];
                    
                    perimeter_force[current_node_global_ci][dim] +=
                            perimeter_force_constant*pow(current_cell_perimeter, 1)*
                            (tension_unit_normals_clockwise[dim]+
                            tension_unit_normals_anti_clockwise[dim]);
                    
                    tension_force[current_node_global_ci][dim] +=
                            clockwise_edge_tension_factor*edge_lengths[current_node_local]*
                            tension_unit_normals_clockwise[dim] +
                            anti_clockwise_edge_tension_factor*edge_lengths[anti_clockwise_node_local]*
                            tension_unit_normals_anti_clockwise[dim];
                    
//                 if(current_cell==12||current_cell==13){
//                     printf("%f ",elongation_force[current_node_global_ci][dim]);
//                 }
                    
// 				printf("%f %f %f \n",tension_force[dim],elongation_force[dim],deformation_force[dim]);
// 				printf("%f \n",internal_angle);
                    
                }
//             if(current_cell==12||current_cell==13){
//                 printf("\n");
//             }
            }
        }
    }
    
    
/////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// boundary force /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
    
    unsigned no_boundary_nodes = mxGetN(prhs[2]);
    
    double target_internal_angle = pi-(2*pi/no_boundary_nodes);
    
    for(unsigned current_node_local=0;current_node_local<no_boundary_nodes;current_node_local++){
        
// 		printf("hello 2 \n");
        
        bool concave_logical = false;
        
        unsigned anti_clockwise_node_local = (current_node_local+no_boundary_nodes-1)%no_boundary_nodes;
        unsigned clockwise_node_local = (current_node_local+1)%no_boundary_nodes;
        
        unsigned current_node_global_mi = (unsigned)boundary_element[current_node_local];
        unsigned anti_clockwise_node_global_mi = (unsigned)boundary_element[anti_clockwise_node_local];
        unsigned clockwise_node_global_mi = (unsigned)boundary_element[clockwise_node_local];
        
        unsigned current_node_global_ci = current_node_global_mi-1;
        unsigned anti_clockwise_node_global_ci = anti_clockwise_node_global_mi-1;
        unsigned clockwise_node_global_ci = clockwise_node_global_mi-1;
        
        double current_node_position[2];
        double anti_clockwise_node_position[2];
        double clockwise_node_position[2];
        
        for(unsigned dim=0;dim<2;dim++) {
            
            unsigned matrix_offset = dim*array_sizes;
            
            current_node_position[dim] = initial_node_positions[current_node_global_ci+matrix_offset];
            anti_clockwise_node_position[dim] = initial_node_positions[anti_clockwise_node_global_ci+matrix_offset];
            clockwise_node_position[dim] = initial_node_positions[clockwise_node_global_ci+matrix_offset];
            
        }
        
        double edge_length_clockwise =
                findStraightLineDistanceBetweenTwoNodes(current_node_position, clockwise_node_position);
        
        double edge_length_anti_clockwise =
                findStraightLineDistanceBetweenTwoNodes(current_node_position, anti_clockwise_node_position);
        
        double deviation_from_average_edge_length_clockwise = edge_length_clockwise-mean_edge_length;
        double deviation_from_average_edge_length_anti_clockwise = edge_length_anti_clockwise-mean_edge_length;
        
        double unit_vector_to_anti_clockwise_node[2];
        double unit_vector_to_clockwise_node[2];
        
        findUnitVector(current_node_position, anti_clockwise_node_position,
                unit_vector_to_anti_clockwise_node);
        findUnitVector(current_node_position, clockwise_node_position,
                unit_vector_to_clockwise_node);
        
        if ((unit_vector_to_anti_clockwise_node[0]*unit_vector_to_clockwise_node[1] -
                unit_vector_to_clockwise_node[0]*unit_vector_to_anti_clockwise_node[1]) < 0) {
            
            concave_logical = true;
        }
        
        double dot_product_unit_vectors = 0;
        double sum_unit_vectors[2];
        
        for(unsigned dim=0;dim<2;dim++) {
            
            sum_unit_vectors[dim] = unit_vector_to_anti_clockwise_node[dim] +
                    unit_vector_to_clockwise_node[dim];
            dot_product_unit_vectors += unit_vector_to_anti_clockwise_node[dim]*
                    unit_vector_to_clockwise_node[dim];
            
        }
        
        if(dot_product_unit_vectors<-1){
            dot_product_unit_vectors = -1;
        }
        else if(dot_product_unit_vectors>1){
            dot_product_unit_vectors = 1;
        }
        
        double internal_angle;
        
        if(concave_logical){
            internal_angle = 2*pi-acos(dot_product_unit_vectors);
        }
        else {
            internal_angle = acos(dot_product_unit_vectors);
        }
        
        double deviation_from_average_angle = internal_angle-target_internal_angle;
        
        double pressure_unit_normals[2];
        
        if(fabs(sum_unit_vectors[0])<0.000000001 && fabs(sum_unit_vectors[1])<0.000000001){
            
            pressure_unit_normals[0] = -unit_vector_to_anti_clockwise_node[1];
            pressure_unit_normals[1] = unit_vector_to_anti_clockwise_node[0];
            
        }
        
        else{
            
            double magnitude_resolved_unit_vectors = findNorm(sum_unit_vectors);
            
            if(concave_logical){
                for(unsigned dim=0;dim<2;dim++) {
                    pressure_unit_normals[dim] =
                            sum_unit_vectors[dim]/magnitude_resolved_unit_vectors;
                }
            }
            else {
                for(unsigned dim=0;dim<2;dim++) {
                    pressure_unit_normals[dim] =
                            -sum_unit_vectors[dim]/magnitude_resolved_unit_vectors;
                }
            }
        }
        
        for(unsigned dim=0;dim<2;dim++) {
            
            boundary_deformation_force[current_node_global_ci][dim] =
                    boundary_deformation_force_constant*fabs(pow(deviation_from_average_angle, 3))/
                    deviation_from_average_angle*pressure_unit_normals[dim];
            
            boundary_edge_force[current_node_global_ci][dim] = boundary_edge_force_constant*(
                    fabs(pow(deviation_from_average_edge_length_clockwise, 3))/
                    deviation_from_average_edge_length_clockwise*
                    unit_vector_to_clockwise_node[dim] +
                    fabs(pow(deviation_from_average_edge_length_anti_clockwise, 3))/
                    deviation_from_average_edge_length_anti_clockwise*
                    unit_vector_to_anti_clockwise_node[dim]);
            
        }
    }
    
/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// apply forces //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
    
    *total_area_force = 0;
    *total_boundary_deformation_force = 0;
    *total_boundary_edge_force = 0;
    *total_deformation_force = 0;
    *total_elongation_force = 0;
    *total_perimeter_force = 0;
    *total_tension_force = 0;
    *total_force = 0;
    
    for(unsigned current_node_global_ci=0;current_node_global_ci<array_sizes;current_node_global_ci++) {
        
        if(cells_per_node[current_node_global_ci]>0){
            
            double vertex_speed[2];
            
            for(unsigned dim=0;dim<2;dim++) {
                
                net_force[current_node_global_ci][dim] +=
                        area_force[current_node_global_ci][dim] +
                        boundary_deformation_force[current_node_global_ci][dim] +
                        boundary_edge_force[current_node_global_ci][dim] +
                        deformation_force[current_node_global_ci][dim] +
                        elongation_force[current_node_global_ci][dim] +
                        perimeter_force[current_node_global_ci][dim] +
                        tension_force[current_node_global_ci][dim];
                
//                 if(current_node_global_ci==248){
//                     printf("%f ", area_force[current_node_global_ci][dim]);
//                 }
                
                unsigned matrix_index = current_node_global_ci+dim*array_sizes;
                
                vertex_speed[dim] = net_force[current_node_global_ci][dim]/viscosity;
                vertex_movements[matrix_index] = vertex_speed[dim]*delta_t;
                
                final_node_positions[matrix_index] =
                        initial_node_positions[matrix_index] + delta_t*vertex_speed[dim];
                
            }
            
//             if(current_node_global_ci==248){
//                 printf("\n");
//             }
            
            double area_force_norm = findNorm(area_force[current_node_global_ci]);
            double boundary_deformation_force_norm = findNorm(boundary_deformation_force[current_node_global_ci]);
            double boundary_edge_force_norm = findNorm(boundary_edge_force[current_node_global_ci]);
            double deformation_force_norm = findNorm(deformation_force[current_node_global_ci]);
            double elongation_force_norm = findNorm(elongation_force[current_node_global_ci]);
            double perimeter_force_norm = findNorm(perimeter_force[current_node_global_ci]);
            double tension_force_norm = findNorm(tension_force[current_node_global_ci]);
            double net_force_norm = findNorm(net_force[current_node_global_ci]);
            
            *total_area_force += area_force_norm;
            *total_boundary_deformation_force += boundary_deformation_force_norm;
            *total_boundary_edge_force += boundary_edge_force_norm;
            *total_deformation_force += deformation_force_norm;
            *total_elongation_force += elongation_force_norm;
            *total_perimeter_force += perimeter_force_norm;
            *total_tension_force += tension_force_norm;
            *total_force += net_force_norm;
            
//             double vertex_speed_magnitude = findNorm(vertex_speed);
            
            
        }
    }
}
