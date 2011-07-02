#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   
   // there must be 3 input arguments
   if (nrhs != 4) {
      mexErrMsgTxt("Exactly 4 input arguments required");
   }
   
   unsigned no_nodes_in_true_solution = mxGetM(prhs[0]);
   
   double* initial_FEM_node_concentration_projection = mxGetPr(prhs[0]);
   double* FEM_edge_nodes_true_solution = mxGetPr(prhs[1]);
   unsigned no_cells = (unsigned)mxGetScalar(prhs[2]);
   unsigned no_nodes_in_test_solution = (unsigned)mxGetScalar(prhs[3]);
   
   plhs[0] = mxCreateDoubleMatrix(no_nodes_in_test_solution, 1, mxREAL);
   double* final_FEM_node_concentration_projection = mxGetPr(plhs[0]);
   
   for(unsigned temp_current_node = 0; temp_current_node<no_nodes_in_test_solution; temp_current_node++){
      final_FEM_node_concentration_projection[temp_current_node] =
              initial_FEM_node_concentration_projection[temp_current_node];
   }
   
   printf("hello \n");
   
   for(unsigned current_node=no_nodes_in_test_solution-no_cells;current_node<no_nodes_in_true_solution-no_cells;current_node++){
      
      unsigned edge_node_1_mi = FEM_edge_nodes_true_solution[current_node];
      unsigned edge_node_2_mi = FEM_edge_nodes_true_solution[current_node+no_nodes_in_true_solution];
      
      unsigned edge_node_1_ci = edge_node_1_mi-1;
      unsigned edge_node_2_ci = edge_node_2_mi-1;
      
      if(FEM_edge_nodes_true_solution[current_node]>0){
         
         final_FEM_node_concentration_projection[current_node] = 0.5*
                 (final_FEM_node_concentration_projection[edge_node_1_ci]+
                 final_FEM_node_concentration_projection[edge_node_2_ci]);
         
      }
      
   }
}

