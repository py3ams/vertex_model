#include<cmath>
#include<cstdlib>
#include "mex.h"
#include "CommonFunctions.hpp"

#define pi 3.1415926535897932384626433832795

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // there must be 4 input arguments
   if (nrhs != 4) {
      mexErrMsgTxt("Exactly 4 input arguments required");
   }
   
   unsigned no_nodes_in_true_solution = mxGetM(prhs[0]);
   
   double* initial_FEM_node_concentration_projection = mxGetPr(prhs[0]);
   double* true_solution_FEM_edge_nodes = mxGetPr(prhs[1]);
   unsigned no_cells = (unsigned)mxGetScalar(prhs[2]);
   unsigned no_nodes_in_test_solution = (unsigned)mxGetScalar(prhs[3]);
   
   plhs[0] = mxCreateDoubleMatrix(no_nodes_in_true_solution, 1, mxREAL);
   double* final_FEM_node_concentration_projection = mxGetPr(plhs[0]);
   
   for(unsigned temp_current_node=0;temp_current_node<no_nodes_in_true_solution;temp_current_node++){
      final_FEM_node_concentration_projection[temp_current_node] =
              initial_FEM_node_concentration_projection[temp_current_node];
   }
   
   for(unsigned current_node=no_nodes_in_test_solution-no_cells;current_node<(no_nodes_in_true_solution-no_cells);current_node++){
      
//       printf("%u \n",current_node);
      
      unsigned edge_node_1_mi = (unsigned)true_solution_FEM_edge_nodes[current_node];
      
      if(edge_node_1_mi>0){

         unsigned edge_node_2_mi = (unsigned)true_solution_FEM_edge_nodes[current_node+no_nodes_in_true_solution];

         unsigned edge_node_1_ci = edge_node_1_mi-1;
         unsigned edge_node_2_ci = edge_node_2_mi-1;

//          printf("%u %u \n",edge_node_1_ci,edge_node_2_ci);
      
         if(true_solution_FEM_edge_nodes[current_node]>0){

            final_FEM_node_concentration_projection[current_node] = 0.5*
                    (final_FEM_node_concentration_projection[edge_node_1_ci]+
                    final_FEM_node_concentration_projection[edge_node_2_ci]);

         }
      }
   }
}
