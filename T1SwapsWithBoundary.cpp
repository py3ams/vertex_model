#include "mex.h"
#include<cmath>
#include<vector>
#include<cassert>
#include "CommonFunctions.hpp"

#define pi 3.1415926535897932384626433832795

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   
   if (nrhs != 15) {
      mexErrMsgTxt("15 input arguments required");
   }
   
   unsigned no_cells = mxGetM(prhs[0]);
   unsigned array_sizes = mxGetM(prhs[1]);
// don't delete again - no_cell_elements is not neccessarily the same as no_cells if FEM_solve_logical is false
   unsigned no_cell_elements = mxGetM(prhs[2]);
   unsigned width_vertex_cells = mxGetN(prhs[3]);
   unsigned no_FEM_elements = mxGetM(prhs[6]);
   unsigned no_FEM_nodes = mxGetM(prhs[8]);
   
   // there are at most 4 edges that can be edited during a single swap, and each could
   // require two changes in the refined_edge_matrix (actually 4 if we consider that (i,j)
   // and (j,i) must be changed for each alteration, but this is taken care of back in
   // matlab. the two changes come from setting one entry to 0 and one entry to 1
   unsigned max_no_refined_edge_matrix_edits = 8;
   
   double* initial_vertex_positions = mxGetPr(prhs[1]);
   double* initial_vertex_cells = mxGetPr(prhs[3]);
   double* initial_cells_per_vertex = mxGetPr(prhs[4]);
   double* initial_Dpp = mxGetPr(prhs[5]);
   double* initial_FEM_elements = mxGetPr(prhs[6]);
   bool* FEM_solve_logical = mxGetLogicals(prhs[7]);
   double* initial_previous_FEM_node_positions = mxGetPr(prhs[8]);
   double protection_time = mxGetScalar(prhs[9]);
   double T1_probability = mxGetScalar(prhs[10]);
   double threshold_T1_swaps = mxGetScalar(prhs[11]);
   double time = mxGetScalar(prhs[12]);
   double* initial_time_vertices_created = mxGetPr(prhs[13]);
   double* initial_FEM_nodes_edge = mxGetPr(prhs[14]);
   
   plhs[0] = mxCreateCellMatrix(no_cells, 1);
   plhs[1] = mxCreateDoubleMatrix(array_sizes, 2, mxREAL);
   plhs[2] = mxCreateCellMatrix(no_cell_elements, 1);
   plhs[3] = mxCreateDoubleMatrix(array_sizes, width_vertex_cells, mxREAL);
   plhs[4] = mxCreateDoubleMatrix(array_sizes, 1, mxREAL);
   plhs[5] = mxCreateDoubleMatrix(no_FEM_nodes, 1, mxREAL);
   plhs[6] = mxCreateDoubleMatrix(no_FEM_elements, 3, mxREAL);
   plhs[7] = mxCreateDoubleMatrix(1, 1, mxREAL);
   plhs[8] = mxCreateDoubleMatrix(no_FEM_nodes, 2, mxREAL);
   plhs[9] = mxCreateDoubleMatrix(array_sizes, 1, mxREAL);
   plhs[10] = mxCreateDoubleMatrix(no_FEM_nodes, 2, mxREAL);
   plhs[11] = mxCreateDoubleMatrix(max_no_refined_edge_matrix_edits, 3, mxREAL);
   
   double* final_vertex_positions = mxGetPr(plhs[1]);
   double* final_vertex_cells = mxGetPr(plhs[3]);
   double* final_cells_per_vertex = mxGetPr(plhs[4]);
   double* final_Dpp = mxGetPr(plhs[5]);
   double* final_FEM_elements = mxGetPr(plhs[6]);
   double* no_T1_swaps_this_iteration = mxGetPr(plhs[7]);
   double* final_previous_FEM_node_positions = mxGetPr(plhs[8]);
   double* final_time_vertices_created = mxGetPr(plhs[9]);
   double* final_FEM_nodes_edge = mxGetPr(plhs[10]);
   double* refined_edge_matrix_edits = mxGetPr(plhs[11]);
   
   *no_T1_swaps_this_iteration = 0;
   
//     printf("hello \n");
   
// set output cells equal to input cells and same for cell elements
   for(unsigned current_cell=0; current_cell<no_cells; current_cell++) {
      
      mxArray* mx_initial_cell_vertices = mxGetCell(prhs[0], current_cell);
      double* initial_cell_vertices = mxGetPr(mx_initial_cell_vertices);
      
      unsigned no_cell_vertices = mxGetN(mx_initial_cell_vertices);
      
      mxArray* mx_final_cell_vertices = mxCreateDoubleMatrix(1, no_cell_vertices, mxREAL);
      double* final_cell_vertices = mxGetPr(mx_final_cell_vertices);
      
      for(unsigned i=0; i<no_cell_vertices;i++){
         final_cell_vertices[i] = initial_cell_vertices[i];
      }
      
      mxSetCell(plhs[0], current_cell, mx_final_cell_vertices);
      
   }
   
//      printf("hello \n");
//
//      printf("%u \n",no_cell_elements);
   
   for(unsigned current_cell_element=0; current_cell_element<no_cell_elements; current_cell_element++) {
      
      mxArray* mx_initial_cell_elements = mxGetCell(prhs[2], current_cell_element);
      double* initial_cell_elements = mxGetPr(mx_initial_cell_elements);
      
      unsigned no_elements = mxGetN(mx_initial_cell_elements);
      
      mxArray* mx_final_cell_elements = mxCreateDoubleMatrix(1, no_elements, mxREAL);
      double* final_cell_elements = mxGetPr(mx_final_cell_elements);
      
      for(unsigned i=0; i<no_elements;i++){
         final_cell_elements[i] = initial_cell_elements[i];
      }
      
      mxSetCell(plhs[2], current_cell_element, mx_final_cell_elements);
      
   }
   
//       printf("hello \n");
   
   // set all other output matrices equal to input matrices
   for(unsigned current_vertex=0;current_vertex<array_sizes;current_vertex++){
      final_cells_per_vertex[current_vertex] = initial_cells_per_vertex[current_vertex];
      final_time_vertices_created[current_vertex] = initial_time_vertices_created[current_vertex];
      for(unsigned dim=0;dim<2;dim++){
         final_vertex_positions[current_vertex+array_sizes*dim] = initial_vertex_positions[current_vertex+array_sizes*dim];
      }
      for(unsigned i=0;i<width_vertex_cells;i++){
         final_vertex_cells[current_vertex+i*array_sizes] = initial_vertex_cells[current_vertex+i*array_sizes];
      }
   }
   
   for(unsigned current_element=0;current_element<no_FEM_elements;current_element++){
      for(unsigned i=0;i<3;i++){
         final_FEM_elements[current_element+no_FEM_elements*i] =
                 initial_FEM_elements[current_element+no_FEM_elements*i];
      }
   }
   
   for(unsigned current_FEM_node=0;current_FEM_node<no_FEM_nodes;current_FEM_node++){
      final_Dpp[current_FEM_node] = initial_Dpp[current_FEM_node];
      for(unsigned dim=0;dim<2;dim++){
         final_previous_FEM_node_positions[current_FEM_node+no_FEM_nodes*dim] =
                 initial_previous_FEM_node_positions[current_FEM_node+no_FEM_nodes*dim];
         final_FEM_nodes_edge[current_FEM_node+no_FEM_nodes*dim] = initial_FEM_nodes_edge[current_FEM_node+no_FEM_nodes*dim];
      }
   }
   
   bool swap_logical = false;
   unsigned no_refined_edge_matrix_edits = 0;
   
//     printf("hello \n");
   
   
   
   // loop over all cells
   for(unsigned current_cell_ci=0; current_cell_ci<no_cells; current_cell_ci++) {
      
// 		printf("start of main loop \n");
      
      unsigned current_cell_mi = current_cell_ci+1;
      
      mxArray* mx_cell_vertices = mxGetCell(plhs[0], current_cell_ci);
      double* cell_vertices = mxGetPr(mx_cell_vertices);
      
      unsigned no_cell_vertices = mxGetN(mx_cell_vertices);
      
      // only proceed if there are more than 3 cell nodes. this loop closes at the end
      // of the program
      if(no_cell_vertices > 3){
         
         // loop over the nodes of the current cell. this loop also goes to the end of the program
         for(unsigned current_vertex_local=0; current_vertex_local<no_cell_vertices; current_vertex_local++) {
            
            // all these variables need to be available within the if(swap_logical) statement later
            unsigned current_vertex_global_mi, current_vertex_global_ci, clockwise_vertex_global_mi, clockwise_vertex_global_ci;
            unsigned no_cell_vertices_cell_with_same_edge;
            unsigned cell_with_same_edge_ci;
            int cell_with_same_edge_mi;
            double current_vertex_position[2], clockwise_vertex_position[2];
            bool cell_with_same_edge_found;
            double* cell_vertices_cell_with_same_edge;
            
            double rand_number = rand()/((double)RAND_MAX);
            swap_logical = false;
            
            // in the next few statements we will decide whether to set swap_logical
            // to true
            
            // only proceed with p(T1_probability)
            if(rand_number < T1_probability){
               
               // find current node and clockwise node in matlab and c indicies
               current_vertex_global_mi = (unsigned)cell_vertices[current_vertex_local];
               current_vertex_global_ci = current_vertex_global_mi-1;
               
               unsigned clockwise_vertex_local = (current_vertex_local+1)%no_cell_vertices;
               
               clockwise_vertex_global_mi = (unsigned)cell_vertices[clockwise_vertex_local];
               clockwise_vertex_global_ci = clockwise_vertex_global_mi - 1;
               
               // extract the positions of the current node and clockwise node from the node positions matrix
               for(unsigned dim=0;dim<2;dim++) {
                  
                  current_vertex_position[dim] = final_vertex_positions[current_vertex_global_ci + dim*array_sizes];
                  clockwise_vertex_position[dim] = final_vertex_positions[clockwise_vertex_global_ci + dim*array_sizes];
               }
               
               // find the current edge length
               double current_edge_length = findStraightLineDistanceBetweenTwoNodes(current_vertex_position,
                       clockwise_vertex_position);
               
               // only proceed if edge length is less than threshold, and nodes have not been created very
               // recently
               if(current_edge_length < threshold_T1_swaps &&
                       (time-final_time_vertices_created[current_vertex_global_ci])>protection_time &&
                       (time-final_time_vertices_created[clockwise_vertex_global_ci])>protection_time){
                  
                  // a swap is now going to happen, unless we find that the cell_with_same_edge exists
                  // and has only three nodes, in which case we will reset swap_logical to false.
                  swap_logical = true;
                  
                  // look for a cell that shares the edge - we need this for a T1 swap
                  cell_with_same_edge_found = false;
                  cell_with_same_edge_mi = -1;
                  
                  unsigned counter_1 = 0;
                  // loop over all cells containing current_vertex_global
                  for(unsigned i=0;i<final_cells_per_vertex[current_vertex_global_ci];i++){
                     
                     unsigned temp_cell_1_mi = (unsigned)final_vertex_cells[current_vertex_global_ci + i*array_sizes];
                     
//                      printf("%u \n",temp_cell_1_mi);
                     
                     // if not current_cell, loop over all cells containing clockwise_vertex_global
                     if(temp_cell_1_mi != current_cell_mi && !cell_with_same_edge_found){
                        for(unsigned j=0;j<final_cells_per_vertex[clockwise_vertex_global_ci];j++){
                           
                           unsigned temp_cell_2_mi = (unsigned)final_vertex_cells[clockwise_vertex_global_ci + j*array_sizes];
                           
                           // if cell containing clockwise_vertex_global is also a cell containing current_vertex_global,
                           // but not current_cell, then it is cell_with_same_edge;
                           if(temp_cell_2_mi == temp_cell_1_mi){
                              cell_with_same_edge_mi = temp_cell_1_mi;
                              cell_with_same_edge_found = true;
                              break;
                           }
                        }
                     }
                  }
                  
                  if(cell_with_same_edge_found){
                     
//                      printf("%u \n",cell_with_same_edge_mi);
                     
                     cell_with_same_edge_ci = cell_with_same_edge_mi-1;
                     
                     mxArray* mx_cell_vertices_cell_with_same_edge = mxGetCell(plhs[0], cell_with_same_edge_ci);
                     cell_vertices_cell_with_same_edge = mxGetPr(mx_cell_vertices_cell_with_same_edge);
                     no_cell_vertices_cell_with_same_edge = mxGetN(mx_cell_vertices_cell_with_same_edge);
                     
                     // cancel swap if cell with same edge has 3 nodes
                     if(no_cell_vertices_cell_with_same_edge == 3){
                        
                        swap_logical = false;
                        
                     }
                  }
               }
            }
            
            // definitely doing a swap!
            if(swap_logical){
               
//                printf("%u \n", current_vertex_global_mi);
//                printf("%u \n", clockwise_vertex_global_mi);
               
//                printf("%u \n",current_cell_mi);
//                printf("%u \n",cell_with_same_edge_mi);
               
//                								printf("Doing a swap \n");
               (*no_T1_swaps_this_iteration)++;
               
               double edge_midpoint[2];
               for(unsigned dim=0;dim<2;dim++) {
                  edge_midpoint[dim] = (current_vertex_position[dim] +
                          clockwise_vertex_position[dim])/2;
               }
               
               double unit_vector_between_vertices[2];
               findUnitVector(current_vertex_position, clockwise_vertex_position, unit_vector_between_vertices);
               
               double new_vertex_1_direction[2], new_vertex_2_direction[2];
               
               new_vertex_1_direction[0] = unit_vector_between_vertices[1];
               new_vertex_1_direction[1] = -unit_vector_between_vertices[0];
               
               double new_vertex_1_position[2], new_vertex_2_position[2];
               
               for(unsigned dim=0;dim<2;dim++){
                  
                  new_vertex_1_position[dim] = edge_midpoint[dim] +
                          0.55*threshold_T1_swaps*new_vertex_1_direction[dim];
                  
               }
               
               // find two unused nodes to put new nodes in
               unsigned new_vertex_1_ci = 0;
               while(final_cells_per_vertex[new_vertex_1_ci]>0){
                  new_vertex_1_ci++;
               }
               unsigned new_vertex_1_mi = new_vertex_1_ci+1;
               
               assert(new_vertex_1_mi <= array_sizes);
               
               for(unsigned dim=0;dim<2;dim++){
                  final_vertex_positions[new_vertex_1_ci+dim*array_sizes] =
                          new_vertex_1_position[dim];
               }
               
               final_time_vertices_created[new_vertex_1_ci] = time;
               
               unsigned new_vertex_2_ci, new_vertex_2_mi;
               
               if(cell_with_same_edge_found){
                  
                  new_vertex_2_direction[0] = -unit_vector_between_vertices[1];
                  new_vertex_2_direction[1] = unit_vector_between_vertices[0];
                  
                  for(unsigned dim=0;dim<2;dim++){
                     
                     new_vertex_2_position[dim] = edge_midpoint[dim] +
                             0.55*threshold_T1_swaps*new_vertex_2_direction[dim];
                  }
                  
                  new_vertex_2_ci = new_vertex_1_ci+1;
                  while(final_cells_per_vertex[new_vertex_2_ci]>0){
                     new_vertex_2_ci++;
                  }
                  new_vertex_2_mi = new_vertex_2_ci+1;
                  
                  assert(new_vertex_2_mi <= array_sizes);
                  
                  // store new node positions
                  for(unsigned dim=0;dim<2;dim++){
                     
                     final_vertex_positions[new_vertex_2_ci+dim*array_sizes] =
                             new_vertex_2_position[dim];
                  }
                  
                  final_time_vertices_created[new_vertex_2_ci] = time;
                  
               }
// 								printf("found and stored new node positions \n");
               
               // set time new nodes created to current time
               
               // edit current cell vertices to contain new vertex 1 instead of previous two vertices
               mxArray* mx_cell_vertices_edited = mxCreateDoubleMatrix(1, no_cell_vertices-1, mxREAL);
               double* cell_vertices_edited = mxGetPr(mx_cell_vertices_edited);
               
               int temp_counter = -1;
               for(unsigned temp_current_vertex_local=0;temp_current_vertex_local<no_cell_vertices;temp_current_vertex_local++){
                  
                  unsigned temp_current_vertex_global_mi = (unsigned)cell_vertices[temp_current_vertex_local];
                  
                  if(temp_current_vertex_global_mi==current_vertex_global_mi){
                     temp_counter++;
                     cell_vertices_edited[temp_counter]=new_vertex_1_mi;
                  }
                  else if(temp_current_vertex_global_mi==clockwise_vertex_global_mi){}
                  else{
                     temp_counter++;
                     cell_vertices_edited[temp_counter]=temp_current_vertex_global_mi;
                  }
               }
               mxSetCell(plhs[0], current_cell_ci, mx_cell_vertices_edited);
               
// 								printf("edited cell_vertices \n");
               
               final_cells_per_vertex[new_vertex_1_ci] = 1;
               final_vertex_cells[new_vertex_1_ci] = current_cell_mi;
               
               if(cell_with_same_edge_found){
                  
                  // edit cell with same edge to contain new vertex 2 instead of previous two vertices
                  mxArray* mx_cell_vertices_cell_with_same_edge_edited = mxCreateDoubleMatrix(1, no_cell_vertices_cell_with_same_edge-1, mxREAL);
                  double* cell_vertices_cell_with_same_edge_edited = mxGetPr(mx_cell_vertices_cell_with_same_edge_edited);
                  
//                   unsigned no_elements = mxGetN(mx_cell_vertices_cell_with_same_edge_edited);
//                   printf("%u \n",no_elements);
                  
                  temp_counter = -1;
                  for(unsigned temp_current_vertex_local=0;temp_current_vertex_local<no_cell_vertices_cell_with_same_edge;temp_current_vertex_local++){
                     
                     unsigned temp_current_vertex_global_mi = (unsigned)cell_vertices_cell_with_same_edge[temp_current_vertex_local];
                     
                     if(temp_current_vertex_global_mi==current_vertex_global_mi){}
                     else if(temp_current_vertex_global_mi==clockwise_vertex_global_mi){
                        temp_counter++;
                        cell_vertices_cell_with_same_edge_edited[temp_counter]=new_vertex_2_mi;
                     }
                     else{
                        temp_counter++;
                        cell_vertices_cell_with_same_edge_edited[temp_counter]=temp_current_vertex_global_mi;
                     }
                  }

//                   unsigned no_elements = mxGetN(mx_cell_vertices_cell_with_same_edge_edited);
//                   printf("%u \n",no_elements);
                  
//                   printf("%u \n",no_cell_vertices_cell_with_same_edge);
//                   printf("%u \n",temp_counter);
                  mxSetCell(plhs[0], cell_with_same_edge_ci, mx_cell_vertices_cell_with_same_edge_edited);
                  
                  
// 								printf("edited cell_vertices_cell_with_same_edge \n");
                  
                  final_cells_per_vertex[new_vertex_2_ci] = 1;
                  final_vertex_cells[new_vertex_2_ci] = cell_with_same_edge_mi;
                  
//                   printf("%u \n",cell_with_same_edge_mi);
                  
               }
//                else{printf("%f \n",time);}
               
               unsigned stretch_cell_1_global_ci, stretch_cell_1_global_mi;
               bool stretch_cell_1_exists = false;

               if((cell_with_same_edge_found&&final_cells_per_vertex[current_vertex_global_ci]==3)||
                       (!cell_with_same_edge_found&&final_cells_per_vertex[current_vertex_global_ci]==2)){
                  
                  stretch_cell_1_exists = true;
                  
                  if(cell_with_same_edge_found){
                     
                     // we are just looping over the cells that share one of the vertices and finding the one that is not
                     // the current cell or the cell with same edge - this is stretch_cell_1_global
                     for(unsigned i=0;i<3;i++){
                        
                        unsigned temp_cell_1_mi = (unsigned)final_vertex_cells[current_vertex_global_ci + i*array_sizes];
                        
                        if(temp_cell_1_mi!=current_cell_mi && temp_cell_1_mi!=cell_with_same_edge_mi){
                           
                           stretch_cell_1_global_mi = temp_cell_1_mi;
                           
                           stretch_cell_1_global_ci = stretch_cell_1_global_mi-1;
                           
                           // edit stretch cell 1 by adding in two new vertices in place of current node global
                           mxArray* mx_cell_vertices_stretch_cell_1 = mxGetCell(plhs[0], stretch_cell_1_global_ci);
                           double* cell_vertices_stretch_cell_1 = mxGetPr(mx_cell_vertices_stretch_cell_1);
                           unsigned no_cell_vertices_stretch_cell_1 = mxGetN(mx_cell_vertices_stretch_cell_1);
                           
                           mxArray* mx_cell_vertices_stretch_cell_1_edited = mxCreateDoubleMatrix(1, no_cell_vertices_stretch_cell_1+1, mxREAL);
                           double* cell_vertices_stretch_cell_1_edited = mxGetPr(mx_cell_vertices_stretch_cell_1_edited);
                           
                           temp_counter = -1;
                           for(unsigned temp_current_vertex_local=0;temp_current_vertex_local<no_cell_vertices_stretch_cell_1;temp_current_vertex_local++){
                              
                              temp_counter++;
                              unsigned temp_current_vertex_global_mi = (unsigned)cell_vertices_stretch_cell_1[temp_current_vertex_local];
                              
                              if(temp_current_vertex_global_mi==current_vertex_global_mi){
                                 cell_vertices_stretch_cell_1_edited[temp_counter] = new_vertex_2_mi;
                                 temp_counter++;
                                 cell_vertices_stretch_cell_1_edited[temp_counter] = new_vertex_1_mi;
                              }
                              else{
                                 cell_vertices_stretch_cell_1_edited[temp_counter] = temp_current_vertex_global_mi;
                              }
                           }
                           
                           mxSetCell(plhs[0], stretch_cell_1_global_ci, mx_cell_vertices_stretch_cell_1_edited);
                        }
                     }
                     
// add stretch cell 1 to cells_per_vertex and vertex_cells for new_vertex_1 and new_vertex_2
                     final_cells_per_vertex[new_vertex_1_ci]++;
                     unsigned matrix_index = (unsigned)(new_vertex_1_ci+array_sizes*(final_cells_per_vertex[new_vertex_1_ci]-1));
                     final_vertex_cells[matrix_index] = stretch_cell_1_global_mi;
                     
                     final_cells_per_vertex[new_vertex_2_ci]++;
                     matrix_index = (unsigned)(new_vertex_2_ci+array_sizes*(final_cells_per_vertex[new_vertex_2_ci]-1));
                     final_vertex_cells[matrix_index] = stretch_cell_1_global_mi;
                     
                  }
                  
                  else{
                     
                     for(unsigned i=0;i<2;i++){
                        
                        unsigned temp_cell_1_mi = (unsigned)final_vertex_cells[current_vertex_global_ci + i*array_sizes];
                        
                        if(temp_cell_1_mi!=current_cell_mi){
                           
                           stretch_cell_1_global_mi = temp_cell_1_mi;
                           
                           stretch_cell_1_global_ci = stretch_cell_1_global_mi-1;
                           
                            // edit stretch cell 1 by adding in two new vertices in place of current node global
                           mxArray* mx_cell_vertices_stretch_cell_1 = mxGetCell(plhs[0], stretch_cell_1_global_ci);
                           double* cell_vertices_stretch_cell_1 = mxGetPr(mx_cell_vertices_stretch_cell_1);
                           unsigned no_cell_vertices_stretch_cell_1 = mxGetN(mx_cell_vertices_stretch_cell_1);
                           
                           mxArray* mx_cell_vertices_stretch_cell_1_edited = mxCreateDoubleMatrix(1, no_cell_vertices_stretch_cell_1, mxREAL);
                           double* cell_vertices_stretch_cell_1_edited = mxGetPr(mx_cell_vertices_stretch_cell_1_edited);
                           
                           for(unsigned temp_current_vertex_local=0;temp_current_vertex_local<no_cell_vertices_stretch_cell_1;temp_current_vertex_local++){
                              
                              unsigned temp_current_vertex_global_mi = (unsigned)cell_vertices_stretch_cell_1[temp_current_vertex_local];
                              
                              if(temp_current_vertex_global_mi==current_vertex_global_mi){
                                 cell_vertices_stretch_cell_1_edited[temp_current_vertex_local] = new_vertex_1_mi;
                              }
                              else{
                                 cell_vertices_stretch_cell_1_edited[temp_current_vertex_local] = temp_current_vertex_global_mi;
                              }
                           }
                           
                           mxSetCell(plhs[0], stretch_cell_1_global_ci, mx_cell_vertices_stretch_cell_1_edited);
                        }
                        //once we have found the stretch_cell we don't need to keep going
                        break;
                     }
                     // add stretch cell 1 to cells_per_vertex and vertex_cells for new_vertex_1 and new_vertex_2
                     final_cells_per_vertex[new_vertex_1_ci]++;
                     unsigned matrix_index = (unsigned)(new_vertex_1_ci+array_sizes*(final_cells_per_vertex[new_vertex_1_ci]-1));
                     final_vertex_cells[matrix_index] = stretch_cell_1_global_mi;
                     
                     
                     
                  }
                  
               }
               
               unsigned stretch_cell_2_global_ci, stretch_cell_2_global_mi;
               bool stretch_cell_2_exists = false;
               
               if((cell_with_same_edge_found&&final_cells_per_vertex[clockwise_vertex_global_ci]==3)||
                       (!cell_with_same_edge_found&&final_cells_per_vertex[clockwise_vertex_global_ci]==2)){
                  
                  stretch_cell_2_exists = true;
                  
                  if(cell_with_same_edge_found){
                     
                     for(unsigned i=0;i<3;i++){
                        
                        unsigned temp_cell_2_mi = (unsigned)final_vertex_cells[clockwise_vertex_global_ci + i*array_sizes];
                        
                        if(temp_cell_2_mi!=current_cell_mi && temp_cell_2_mi!=cell_with_same_edge_mi){
                           
                           stretch_cell_2_global_mi = temp_cell_2_mi;
                           stretch_cell_2_global_ci = stretch_cell_2_global_mi-1;
                           
                           // edit stretch cell 2 by adding in two new vertices in place of clockwise node global
                           mxArray* mx_cell_vertices_stretch_cell_2 = mxGetCell(plhs[0], stretch_cell_2_global_ci);
                           double* cell_vertices_stretch_cell_2 = mxGetPr(mx_cell_vertices_stretch_cell_2);
                           unsigned no_cell_vertices_stretch_cell_2 = mxGetN(mx_cell_vertices_stretch_cell_2);
                           
                           mxArray* mx_cell_vertices_stretch_cell_2_edited = mxCreateDoubleMatrix(1, no_cell_vertices_stretch_cell_2+1, mxREAL);
                           double* cell_vertices_stretch_cell_2_edited = mxGetPr(mx_cell_vertices_stretch_cell_2_edited);
                           
                           temp_counter = -1;
                           for(unsigned temp_current_vertex_local=0;temp_current_vertex_local<no_cell_vertices_stretch_cell_2;temp_current_vertex_local++){
                              
                              temp_counter++;
                              unsigned temp_current_vertex_global_mi = (unsigned)cell_vertices_stretch_cell_2[temp_current_vertex_local];
                              
//                                                 printf("%u \n", temp_current_vertex_global_mi);
                              
                              if(temp_current_vertex_global_mi==clockwise_vertex_global_mi){
                                 cell_vertices_stretch_cell_2_edited[temp_counter] = new_vertex_1_mi;
                                 temp_counter++;
                                 cell_vertices_stretch_cell_2_edited[temp_counter] = new_vertex_2_mi;
                              }
                              else{
                                 cell_vertices_stretch_cell_2_edited[temp_counter] = temp_current_vertex_global_mi;
                              }
                           }
                           mxSetCell(plhs[0], stretch_cell_2_global_ci, mx_cell_vertices_stretch_cell_2_edited);
                        }
                     }
                     
                     // add stretch cell 2 to cells_per_vertex and vertex_cells for new_vertex_1 and new_vertex_2
                     final_cells_per_vertex[new_vertex_1_ci]++;
                     unsigned matrix_index = (unsigned)(new_vertex_1_ci+array_sizes*(final_cells_per_vertex[new_vertex_1_ci]-1));
                     final_vertex_cells[matrix_index] = stretch_cell_2_global_mi;
                     
                     final_cells_per_vertex[new_vertex_2_ci]++;
                     matrix_index = (unsigned)(new_vertex_2_ci+array_sizes*(final_cells_per_vertex[new_vertex_2_ci]-1));
                     final_vertex_cells[matrix_index] = stretch_cell_2_global_mi;
                  }
                  
                  else{
                     
                     for(unsigned i=0;i<2;i++){
                        
                        unsigned temp_cell_2_mi = (unsigned)final_vertex_cells[clockwise_vertex_global_ci + i*array_sizes];
                        
                        if(temp_cell_2_mi!=current_cell_mi){
                           
                           stretch_cell_2_global_mi = temp_cell_2_mi;
                           stretch_cell_2_global_ci = stretch_cell_2_global_mi-1;
                           
                           // edit stretch cell 2 by adding in new vertex in place of clockwise node global
                           mxArray* mx_cell_vertices_stretch_cell_2 = mxGetCell(plhs[0], stretch_cell_2_global_ci);
                           double* cell_vertices_stretch_cell_2 = mxGetPr(mx_cell_vertices_stretch_cell_2);
                           unsigned no_cell_vertices_stretch_cell_2 = mxGetN(mx_cell_vertices_stretch_cell_2);
                           
                           mxArray* mx_cell_vertices_stretch_cell_2_edited = mxCreateDoubleMatrix(1, no_cell_vertices_stretch_cell_2, mxREAL);
                           double* cell_vertices_stretch_cell_2_edited = mxGetPr(mx_cell_vertices_stretch_cell_2_edited);
                           
                           for(unsigned temp_current_vertex_local=0;temp_current_vertex_local<no_cell_vertices_stretch_cell_2;temp_current_vertex_local++){
                              
                              unsigned temp_current_vertex_global_mi = (unsigned)cell_vertices_stretch_cell_2[temp_current_vertex_local];
                              
//                                                 printf("%u \n", temp_current_vertex_global_mi);
                              
                              if(temp_current_vertex_global_mi==clockwise_vertex_global_mi){
                                 cell_vertices_stretch_cell_2_edited[temp_current_vertex_local] = new_vertex_1_mi;
                              }
                              else{
                                 cell_vertices_stretch_cell_2_edited[temp_current_vertex_local] = temp_current_vertex_global_mi;
                              }
                           }
                           mxSetCell(plhs[0], stretch_cell_2_global_ci, mx_cell_vertices_stretch_cell_2_edited);
                        }
                        // we have found stretch cell 2
                        break;
                     }
                     
                     // add stretch cell 2 to cells_per_vertex and vertex_cells for new_vertex_1 and new_vertex_2
                     final_cells_per_vertex[new_vertex_1_ci]++;
                     unsigned matrix_index = (unsigned)(new_vertex_1_ci+array_sizes*(final_cells_per_vertex[new_vertex_1_ci]-1));
                     final_vertex_cells[matrix_index] = stretch_cell_2_global_mi;
                     
                  }
               }
            
              
               // reset cells per vertex and vertex cells for current_vertex_global and clockwise_vertex_global.
               // they should no longer be part of any cells, as we have replaced them with the new vertices.
               // also set their vertex position to 0. need to be careful we don't try and use these positions anywhere
               // in the FEM stuff. we haven't reset the correspoding FEM_node_positions so these should be used instead.
               final_cells_per_vertex[current_vertex_global_ci] = 0;
               final_cells_per_vertex[clockwise_vertex_global_ci] = 0;
               
               for(unsigned i=0;i<width_vertex_cells;i++){
                  final_vertex_cells[current_vertex_global_ci+i*array_sizes] = 0;
                  final_vertex_cells[clockwise_vertex_global_ci+i*array_sizes] = 0;
                  
               }
               
               for(unsigned dim=0;dim<2;dim++){
                  final_vertex_positions[current_vertex_global_ci+dim*array_sizes] = 0;
                  final_vertex_positions[clockwise_vertex_global_ci+dim*array_sizes] = 0;
               }
               
//                // at this point the cells should be completely edited. we can now do all the FEM stuff
//                if(*FEM_solve_logical){
//
//                   unsigned new_node_1_ci = new_vertex_1_ci;
//                   unsigned new_node_2_ci = new_vertex_2_ci;
//
//                   unsigned new_node_1_mi = new_vertex_1_mi;
//                   unsigned new_node_2_mi = new_vertex_2_mi;
//
//                   for(unsigned dim=0;dim<2;dim++){
//
//                      final_previous_FEM_node_positions[new_node_1_ci+dim*no_FEM_nodes] = new_vertex_1_position[dim];
//                      if(cell_with_same_edge_found){
//                         final_previous_FEM_node_positions[new_node_2_ci+dim*no_FEM_nodes] = new_vertex_2_position[dim];
//                      }
//                   }
//
//                   // find Dpp at edge mid-point as well as centroids of the two cells that share the edge
//                   double Dpp_at_old_edge_mid_point = (final_Dpp[current_vertex_global_ci] +
//                           final_Dpp[clockwise_vertex_global_ci])/2;
//
//                   final_Dpp[new_vertex_1_ci] = Dpp_at_old_edge_mid_point;
//                   if(cell_with_same_edge_found){
//                      final_Dpp[new_vertex_2_ci] = Dpp_at_old_edge_mid_point;
//                   }
//
//                   // there will be 2 reusable elements - one each from the current cell and cell_with_same_edge.
//                   // we can use these for the new elements that will appear in the stretch cells
//                   unsigned reusable_elements_ci[2];
//
//                   // edit the elements of the current cell
//
//                   mxArray* mx_current_cell_elements = mxGetCell(plhs[2], current_cell_ci);
//                   double* current_cell_elements_mi = mxGetPr(mx_current_cell_elements);
//                   unsigned no_elements_current_cell = mxGetN(mx_current_cell_elements);
//
//                   mxArray* mx_current_cell_elements_edited = mxCreateDoubleMatrix(1, no_elements_current_cell-1, mxREAL);
//                   double* current_cell_elements_edited = mxGetPr(mx_current_cell_elements_edited);
//
//                   // loop over the elements associated with the current cell. we will need to edit 3 elements in this cell.
//                   // one will be completely removed, and 2 will have a node changed from the old node to the new node.
//                   temp_counter = -1;
//                   for(unsigned i=0;i<no_elements_current_cell;i++){
//
//                      unsigned current_element_mi = (unsigned)current_cell_elements_mi[i];
//                      unsigned current_element_ci = current_element_mi-1;
//
//                      bool remove_this_element = false;
//
//                      // loop over the nodes of the current element to find if it contains any of the old vertices
//                      for(unsigned j=0; j<3; j++){
//
//                         unsigned temp_current_FEM_node_mi = (unsigned)final_FEM_elements[current_element_ci+j*no_FEM_elements];
//                         unsigned temp_clockwise_FEM_node_mi = (unsigned)final_FEM_elements[current_element_ci+((j+1)%3)*no_FEM_elements];
//
//                         unsigned temp_anti_clockwise_FEM_node_mi = (unsigned)final_FEM_elements[current_element_ci+((j-1)%3)*no_FEM_elements];
//                         unsigned temp_anti_clockwise_FEM_node_ci = temp_anti_clockwise_FEM_node_mi-1;
//
//                         // if current element contains both the old nodes we are going to remove this element and eventually re-use it
//                         // for one of the new elements in the stretch cells
//                         if(temp_current_FEM_node_mi==current_vertex_global_mi &&
//                                 temp_clockwise_FEM_node_mi == clockwise_vertex_global_mi){
//
//                            remove_this_element = true;
//
//                            final_FEM_elements[current_element_ci+j*no_FEM_elements] = 0;
//                            final_FEM_elements[current_element_ci+((j+1)%3)*no_FEM_elements] = 0;
//                            final_FEM_elements[current_element_ci+((j+2)%3)*no_FEM_elements] = 0;
//
//                            reusable_elements_ci[0] = current_element_ci;
//
//                            break;
//                         }
//
//                         // if the element contains instead just one old node we need to replace it with the appropriate new node
//
//                         else if(temp_current_FEM_node_mi==current_vertex_global_mi &&
//                                 temp_clockwise_FEM_node_mi==(no_FEM_nodes-no_cells+current_cell_mi)){
//
//                            final_FEM_elements[current_element_ci+j*no_FEM_elements] = new_node_1_mi;
//
//                            // this if statement is somewhat redundant, but adds a bit of clarity to some complicated code
//                            if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci]>0){
//
//                               // if the third node in this element is an edge node, we must edit the edge it is on to contain new_vertex_1 rather
//                               // than the old current_vertex. as FEM_nodes_edge can contain the edge in either order, we must check if current_vertex
//                               // is in the first position or the second position. this process was originally implemented
//                               // later in the code, during the loops over the stretch cells. in most cases these are equivalent, however in rare cases
//                               // the stretch cells do not exist, so it is important to do it here with the current_cell and cell_with_same_edge,
//                               // which must both exist during a T1 swap.
//
//                               if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci]==current_vertex_global_mi){
//
//                                  final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci] = new_vertex_1_mi;
//
//                                  // refined_edge_matrix is sparse in Matlab, so refined_edge_matrix_edits is an Nx3 matrix with the
//                                  // row and column index of the change then a 0 or 1.
//                                  no_refined_edge_matrix_edits++;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = current_vertex_global_mi;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                          final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
//
//                                  no_refined_edge_matrix_edits++;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_1_mi;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                          final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
//
//                               }
//
//                               // this is the case where current_vertex_global happens to be in the second column of FEM_nodes_edge
//                               // at the temp_anti_clockwise_FEM_node_ci.
//
//                               else if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes]==current_vertex_global_mi){
//
//                                  final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes] = new_vertex_1_mi;
//
//                                  no_refined_edge_matrix_edits++;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = current_vertex_global_mi;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                          final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
//
//                                  no_refined_edge_matrix_edits++;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_1_mi;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                          final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
//
//                               }
//                            }
//
//                            break;
//
//                         }
//
//                         else if(temp_current_FEM_node_mi==(no_FEM_nodes-no_cells+current_cell_mi) &&
//                                 temp_clockwise_FEM_node_mi==clockwise_vertex_global_mi){
//
//                            final_FEM_elements[current_element_ci+((j+1)%3)*no_FEM_elements] = new_node_1_mi;
//
//                            // check if the third node is an edge node
//                            if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci]>0){
//
//                               // check whether clockwise_vertex_global is in the first column of the edge
//                               if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci]==clockwise_vertex_global_mi){
//
//                                  final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci] = new_vertex_1_mi;
//
//                                  no_refined_edge_matrix_edits++;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = clockwise_vertex_global_mi;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                          final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
//
//                                  no_refined_edge_matrix_edits++;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_1_mi;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                          final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
//
//                               }
//
//                               // alternatively check if clockwise_vertex_global is in the second column of the edge. again this line
//                               // could be done away with, as if the node is not in the first column it should be in the second column.
//                               // we have already checked that the node is an edge node.
//                               else if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes]==clockwise_vertex_global_mi){
//
//                                  final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes] = new_vertex_1_mi;
//
//                                  no_refined_edge_matrix_edits++;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = clockwise_vertex_global_mi;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                          final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
//
//                                  no_refined_edge_matrix_edits++;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_1_mi;
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                          final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
//                                  refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
//
//                               }
//                            }
//                         }
//                      }
//
//                      // we don't want to mess around with the indexing of the for loop so we use this counter
//                      // to store new elements when not removing a ndoe
//                      if(!remove_this_element){
//
//                         temp_counter++;
//                         current_cell_elements_edited[temp_counter] = current_element_mi;
//
//                      }
//                   }
//
//                   mxSetCell(plhs[2], current_cell_ci, mx_current_cell_elements_edited);
//
//                   // now, let's do the same with cell_with_same_edge
//
//                   if(cell_with_same_edge_found){
//                      mxArray* mx_cell_with_same_remove_this_elements = mxGetCell(plhs[2], cell_with_same_edge_ci);
//                      double* cell_with_same_remove_this_elements_mi = mxGetPr(mx_cell_with_same_remove_this_elements);
//
//                      unsigned no_elements_cell_with_same_edge = mxGetN(mx_cell_with_same_remove_this_elements);
//
//                      mxArray* mx_cell_with_same_remove_this_elements_edited = mxCreateDoubleMatrix(1, no_elements_cell_with_same_edge-1, mxREAL);
//                      double* cell_with_same_remove_this_elements_edited = mxGetPr(mx_cell_with_same_remove_this_elements_edited);
//
//                      temp_counter = -1;
//                      for(unsigned i=0;i<no_elements_cell_with_same_edge;i++){
//
//                         unsigned current_element_mi = (unsigned)cell_with_same_remove_this_elements_mi[i];
//                         unsigned current_element_ci = current_element_mi-1;
//
//                         bool remove_this_element = false;
//
//                         for(unsigned j=0; j<3; j++){
//
//                            unsigned temp_current_FEM_node_mi =
//                                    (unsigned)final_FEM_elements[current_element_ci+j*no_FEM_elements];
//                            unsigned temp_clockwise_FEM_node_mi =
//                                    (unsigned)final_FEM_elements[current_element_ci+((j+1)%3)*no_FEM_elements];
//
//                            unsigned temp_anti_clockwise_FEM_node_mi =
//                                    (unsigned)final_FEM_elements[current_element_ci+((j-1)%3)*no_FEM_elements];
//                            unsigned temp_anti_clockwise_FEM_node_ci = temp_anti_clockwise_FEM_node_mi-1;
//
//                            if(temp_current_FEM_node_mi==clockwise_vertex_global_mi &&
//                                    temp_clockwise_FEM_node_mi == current_vertex_global_mi){
//
//                               remove_this_element = true;
//
//                               final_FEM_elements[current_element_ci+j*no_FEM_elements] = 0;
//                               final_FEM_elements[current_element_ci+((j+1)%3)*no_FEM_elements] = 0;
//                               final_FEM_elements[current_element_ci+((j+2)%3)*no_FEM_elements] = 0;
//
//                               reusable_elements_ci[1] = current_element_ci;
//
//                               break;
//                            }
//
//                            // may need to play with FEM_nodes_edge and refined_edge_matrix here
//
//                            else if(temp_current_FEM_node_mi==clockwise_vertex_global_mi &&
//                                    temp_clockwise_FEM_node_mi==(no_FEM_nodes-no_cells+cell_with_same_edge_mi)){
//
//                               final_FEM_elements[current_element_ci+j*no_FEM_elements] = new_vertex_2_mi;
//
//                               // see above (comments in current_cell) for what is going on here
//                               if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci]>0){
//
//                                  if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci]==clockwise_vertex_global_mi){
//
//                                     final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci] = new_vertex_2_mi;
//
//                                     no_refined_edge_matrix_edits++;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = clockwise_vertex_global_mi;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                             final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
//
//                                     no_refined_edge_matrix_edits++;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_2_mi;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                             final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
//
//                                  }
//
//                                  else if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes]==clockwise_vertex_global_mi){
//
//                                     final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes] = new_vertex_2_mi;
//
//                                     no_refined_edge_matrix_edits++;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = clockwise_vertex_global_mi;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                             final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
//
//                                     no_refined_edge_matrix_edits++;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_2_mi;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                             final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
//
//                                  }
//                               }
//                               break;
//
//                            }
//
//                            else if(temp_current_FEM_node_mi==(no_FEM_nodes-no_cells+cell_with_same_edge_mi) &&
//                                    temp_clockwise_FEM_node_mi==current_vertex_global_mi){
//
//                               final_FEM_elements[current_element_ci+((j+1)%3)*no_FEM_elements] = new_vertex_2_mi;
//
//                               if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci]>0){
//
//                                  if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci]==current_vertex_global_mi){
//
//                                     final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci] = new_vertex_2_mi;
//
//                                     no_refined_edge_matrix_edits++;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = current_vertex_global_mi;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                             final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
//
//                                     no_refined_edge_matrix_edits++;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_2_mi;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                             final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
//
//                                  }
//
//                                  // see comments above about what is going on here. this is the case where current_vertex_global happens
//                                  // to be in the second column of FEM_nodes_edge at the temp_anti_clockwise_FEM_node_ci.
//                                  else if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes]==current_vertex_global_mi){
//
//                                     final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes] = new_vertex_2_mi;
//
//                                     no_refined_edge_matrix_edits++;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = current_vertex_global_mi;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                             final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
//
//                                     no_refined_edge_matrix_edits++;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_2_mi;
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
//                                             final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
//                                     refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
//
//                                  }
//                               }
//                            }
//                         }
//
//                         if(!remove_this_element){
//
//                            temp_counter++;
//                            cell_with_same_remove_this_elements_edited[temp_counter] = current_element_mi;
//
//                         }
//
//                      }
//
//                      mxSetCell(plhs[2], cell_with_same_edge_ci, mx_cell_with_same_remove_this_elements_edited);
//                   }
//
//                   if(stretch_cell_1_exists){
//
//                      mxArray* mx_stretch_cell_1_elements = mxGetCell(plhs[2], stretch_cell_1_global_ci);
//                      double* stretch_cell_1_elements_mi = mxGetPr(mx_stretch_cell_1_elements);
//                      unsigned no_elements_stretch_cell_1 = mxGetN(mx_stretch_cell_1_elements);
//
//                      if(cell_with_same_edge_found){
//
//                         mxArray* mx_stretch_cell_1_elements_edited = mxCreateDoubleMatrix(1, no_elements_stretch_cell_1+1, mxREAL);
//                         double* stretch_cell_1_elements_edited = mxGetPr(mx_stretch_cell_1_elements_edited);
//
//                         temp_counter = -1;
//
//                         // loop over the elements of stretch_cell_1
//                         for(unsigned i=0;i<no_elements_stretch_cell_1;i++){
//
//                            unsigned current_element_mi = (unsigned)stretch_cell_1_elements_mi[i];
//                            unsigned current_element_ci = current_element_mi-1;
//
//                            bool new_element_now = false;
//
//                            for(unsigned j=0; j<3; j++){
//
//                               unsigned temp_current_FEM_node_mi = (unsigned)final_FEM_elements[current_element_ci+j*no_FEM_elements];
//                               unsigned temp_clockwise_FEM_node_mi = (unsigned)final_FEM_elements[current_element_ci+((j+1)%3)*no_FEM_elements];
//
//                               unsigned temp_anti_clockwise_FEM_node_mi = (unsigned)final_FEM_elements[current_element_ci+((j-1)%3)*no_FEM_elements];
//                               unsigned temp_anti_clockwise_FEM_node_ci = temp_anti_clockwise_FEM_node_mi-1;
//
//                               // we are looking to see if the current element is the one containing current_vertex_global with the central node
//                               // of stretch_cell_1 clockwise to it. in this instance the current_vertex must be replaced by new_vertex_2 (probably
//                               // need to draw out a T1 swap to see why this is so.
//                               if(temp_current_FEM_node_mi==current_vertex_global_mi &&
//                                       temp_clockwise_FEM_node_mi==(no_FEM_nodes-no_cells+stretch_cell_1_global_mi)){
//
//                                  final_FEM_elements[current_element_ci+j*no_FEM_elements] = new_vertex_2_mi;
//
// //                                                     // if the third node in this element is an edge node, we must edit the edge it is on to contain new_vertex_2 rather
// //                                                     // than the old current_vertex. as FEM_nodes_edge can contain the edge in either order, we must check if current_vertex
// //                                                     // is in the first position or the second position. we could first check whether either entry is non-zero, though this
// //                                                     // would just add an unneccesary if statement (though may be clearer what's going on). it should also be noted that
// //                                                     // this process could equally have been carried out when we were editing the elements of current_cell or cell_with_same_edge
// //                                                     // as these edges are shared with those cells.
// //                                                     if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci]==current_vertex_global_mi){
// //
// //                                                         final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci] = new_vertex_2_mi;
// //
// //                                                         // we need to remove the old FEM_edge from the refined_edge_matrix, and add in the new one instead.
// //                                                         // refined_edge_matrix is sparse in Matlab, so refined_edge_matrix_edits is an Nx3 matrix with the
// //                                                         // row and column index of the change then a 0 or 1.
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = current_vertex_global_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
// //
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_2_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
// //
// //                                                     }
// //
// //                                                     // see comments above about what is going on here. this is the case where current_vertex_global happens
// //                                                     // to be in the second column of FEM_nodes_edge at the temp_anti_clockwise_FEM_node_ci.
// //                                                     else if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes]==current_vertex_global_mi){
// //
// //                                                         final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes] = new_vertex_2_mi;
// //
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = current_vertex_global_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
// //
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_2_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
// //
// //                                                     }
//
//                                  new_element_now = true;
//
//                                  break;
//
//                               }
//
//                               else if(temp_current_FEM_node_mi==(no_FEM_nodes-no_cells+stretch_cell_1_global_mi) &&
//                                       temp_clockwise_FEM_node_mi==current_vertex_global_mi){
//
//                                  final_FEM_elements[current_element_ci+((j+1)%3)*no_FEM_elements] = new_vertex_1_mi;
//
// //                                                     // We need to edit final_FEM_nodes_edge to contain the new nodes if the new nodes connect to an edge containting an FEM node.
// //                                                     if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci]==current_vertex_global_mi){
// //
// //                                                         final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci] = new_vertex_1_mi;
// //
// //                                                         // see comments above
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = current_vertex_global_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
// //
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_1_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
// //
// //                                                     }
// //
// //                                                     else if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes]==current_vertex_global_mi){
// //
// //                                                         final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes] = new_vertex_1_mi;
// //
// //                                                         // see comments above
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = current_vertex_global_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
// //
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_1_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
// //
// //                                                     }
//
//                                  break;
//
//                               }
//                            }
//
//                            temp_counter++;
//                            stretch_cell_1_elements_edited[temp_counter] = current_element_mi;
//
//                            if(new_element_now){
//
//                               temp_counter++;
//
//                               unsigned new_element_ci = reusable_elements_ci[0];
//                               unsigned new_element_mi = new_element_ci + 1;
//
//                               final_FEM_elements[new_element_ci] = new_vertex_2_mi;
//                               final_FEM_elements[new_element_ci+no_FEM_elements] = new_vertex_1_mi;
//                               final_FEM_elements[new_element_ci+2*no_FEM_elements] = no_FEM_nodes-no_cells+stretch_cell_1_global_mi;
//
//                               stretch_cell_1_elements_edited[temp_counter] = new_element_mi;
//
//                            }
//                         }
//
//                         mxSetCell(plhs[2], stretch_cell_1_global_ci, mx_stretch_cell_1_elements_edited);
//                      }
//                      else{
//
//                         for(unsigned i=0;i<no_elements_stretch_cell_1;i++){
//
//                            unsigned current_element_mi = (unsigned)stretch_cell_1_elements_mi[i];
//                            unsigned current_element_ci = current_element_mi-1;
//
//                            for(unsigned j=0; j<3; j++){
//
//                               unsigned temp_current_FEM_node_mi = (unsigned)final_FEM_elements[current_element_ci+j*no_FEM_elements];
//
//                               if(temp_current_FEM_node_mi==current_vertex_global_mi){
//
//                                  final_FEM_elements[current_element_ci+j*no_FEM_elements] = new_vertex_1_mi;
//
//                               }
//                            }
//                         }
//                      }
//                   }
//
//
//                   if(stretch_cell_2_exists){
//
//                      if(cell_with_same_edge_found){
//
//                         mxArray* mx_stretch_cell_2_elements = mxGetCell(plhs[2], stretch_cell_2_global_ci);
//                         double* stretch_cell_2_elements_mi = mxGetPr(mx_stretch_cell_2_elements);
//
//                         unsigned no_elements_stretch_cell_2 = mxGetN(mx_stretch_cell_2_elements);
//
//                         if(cell_with_same_edge_found){
//
//                            mxArray* mx_stretch_cell_2_elements_edited = mxCreateDoubleMatrix(1, no_elements_stretch_cell_2+1, mxREAL);
//                            double* stretch_cell_2_elements_edited = mxGetPr(mx_stretch_cell_2_elements_edited);
//
//                            temp_counter = -1;
//                            for(unsigned i=0;i<no_elements_stretch_cell_2;i++){
//
//                               unsigned current_element_mi = (unsigned)stretch_cell_2_elements_mi[i];
//                               unsigned current_element_ci = current_element_mi-1;
//
//                               bool new_element_now = false;
//
//                               for(unsigned j=0; j<3; j++){
//
//                                  unsigned temp_current_FEM_node_mi = (unsigned)final_FEM_elements[current_element_ci+j*no_FEM_elements];
//                                  unsigned temp_clockwise_FEM_node_mi = (unsigned)final_FEM_elements[current_element_ci+((j+1)%3)*no_FEM_elements];
//
//                                  unsigned temp_anti_clockwise_FEM_node_mi = (unsigned)final_FEM_elements[current_element_ci+((j-1)%3)*no_FEM_elements];
//                                  unsigned temp_anti_clockwise_FEM_node_ci = temp_anti_clockwise_FEM_node_mi-1;
//
//                                  if(temp_current_FEM_node_mi==clockwise_vertex_global_mi &&
//                                          temp_clockwise_FEM_node_mi==(no_FEM_nodes-no_cells+stretch_cell_2_global_mi)){
//
//                                     final_FEM_elements[current_element_ci+j*no_FEM_elements] = new_vertex_1_mi;
//
//
// //                                                     // Comments about what is going on here can be found above for stretch_cell_1
// //                                                     if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci]==clockwise_vertex_global_mi){
// //
// //                                                         final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci] = new_vertex_1_mi;
// //
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = clockwise_vertex_global_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
// //
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_1_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
// //
// //                                                     }
// //
// //                                                     else if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes]==clockwise_vertex_global_mi){
// //
// //                                                         final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes] = new_vertex_1_mi;
// //
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = clockwise_vertex_global_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
// //
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_1_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
// //
// //                                                     }
//
//                                     new_element_now = true;
//
//                                     break;
//
//                                  }
//
//                                  else if(temp_current_FEM_node_mi==(no_FEM_nodes-no_cells+stretch_cell_2_global_mi) &&
//                                          temp_clockwise_FEM_node_mi==clockwise_vertex_global_mi){
//
//                                     final_FEM_elements[current_element_ci+((j+1)%3)*no_FEM_elements] = new_vertex_2_mi;
//
// //                                                     // We need to edit final_FEM_nodes_edge to contain the new nodes if the new nodes connect to an edge containting an FEM node.
// //                                                     if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci]==clockwise_vertex_global_mi){
// //
// //                                                         final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci] = new_vertex_2_mi;
// //
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = clockwise_vertex_global_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
// //
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_2_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
// //
// //                                                     }
// //
// //                                                     else if(final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes]==clockwise_vertex_global_mi){
// //
// //                                                         final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci+no_FEM_nodes] = new_vertex_2_mi;
// //
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = clockwise_vertex_global_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 0;
// //
// //                                                         no_refined_edge_matrix_edits++;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1] = new_vertex_2_mi;
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+max_no_refined_edge_matrix_edits] =
// //                                                                 final_FEM_nodes_edge[temp_anti_clockwise_FEM_node_ci];
// //                                                         refined_edge_matrix_edits[no_refined_edge_matrix_edits-1+2*max_no_refined_edge_matrix_edits] = 1;
// //
// //                                                     }
//
//                                     break;
//
//                                  }
//                               }
//
//                               temp_counter++;
//                               stretch_cell_2_elements_edited[temp_counter] = current_element_mi;
//
//                               if(new_element_now){
//
//                                  temp_counter++;
//
//                                  unsigned new_element_ci = reusable_elements_ci[1];
//                                  unsigned new_element_mi = new_element_ci + 1;
//
//                                  final_FEM_elements[new_element_ci] = new_vertex_1_mi;
//                                  final_FEM_elements[new_element_ci+no_FEM_elements] = new_vertex_2_mi;
//                                  final_FEM_elements[new_element_ci+2*no_FEM_elements] = no_FEM_nodes-no_cells+stretch_cell_2_global_mi;
//
//                                  stretch_cell_2_elements_edited[temp_counter] = new_element_mi;
//
//                               }
//                            }
//
//                            mxSetCell(plhs[2], stretch_cell_2_global_ci, mx_stretch_cell_2_elements_edited);
//                         }
//                         else{
//
//                            for(unsigned i=0;i<no_elements_stretch_cell_2;i++){
//
//                               unsigned current_element_mi = (unsigned)stretch_cell_2_elements_mi[i];
//                               unsigned current_element_ci = current_element_mi-1;
//
//                               for(unsigned j=0; j<3; j++){
//
//                                  unsigned temp_current_FEM_node_mi = (unsigned)final_FEM_elements[current_element_ci+j*no_FEM_elements];
//
//                                  if(temp_current_FEM_node_mi==clockwise_vertex_global_mi){
//
//                                     final_FEM_elements[current_element_ci+j*no_FEM_elements] = new_vertex_1_mi;
//
//                                  }
//                               }
//                            }
//
//                         }
//
//                         for(unsigned dim=0;dim<2;dim++){
//
//                            final_previous_FEM_node_positions[new_vertex_1_ci+dim*no_FEM_nodes] = new_vertex_1_position[dim];
//                            if(cell_with_same_edge_found){
//                               final_previous_FEM_node_positions[new_vertex_2_ci+dim*no_FEM_nodes] = new_vertex_2_position[dim];
//                            }
//
//                            final_previous_FEM_node_positions[current_vertex_global_ci+dim*no_FEM_nodes] = 0;
//                            final_previous_FEM_node_positions[clockwise_vertex_global_ci+dim*no_FEM_nodes] = 0;
//
//                         }
//
//                         final_Dpp[current_vertex_global_ci] = 0;
//                         final_Dpp[clockwise_vertex_global_ci] = 0;
//                      }
//                   }
//                }
               break;
            }
         }
      }
      if(swap_logical){break;}
   }
}
