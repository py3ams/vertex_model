#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    double* cell_store = mxGetPr(prhs[0]);
    double* cells_per_node = mxGetPr(prhs[2]);
    
    unsigned array_sizes = mxGetM(prhs[2]);
    
    unsigned no_boundary_nodes = 0;
    
    for(unsigned current_node_ci=0;current_node_ci<array_sizes;current_node_ci++){
        if(cells_per_node[current_node_ci]>0.5 && cells_per_node[current_node_ci]<2.5){
            no_boundary_nodes++;
        }
    }
    
    plhs[0] = mxCreateDoubleMatrix(1, no_boundary_nodes, mxREAL);
    double* boundary_element = mxGetPr(plhs[0]);
    
    unsigned current_boundary_node_ci;
    
    for(unsigned current_node_ci=0;current_node_ci<array_sizes;current_node_ci++){
        if(cells_per_node[current_node_ci]>0.5 && cells_per_node[current_node_ci]<1.5){
            current_boundary_node_ci = current_node_ci;
            break;
        }
    }
    
    unsigned current_boundary_node_mi = current_boundary_node_ci+1;
    boundary_element[0] = current_boundary_node_mi;
    
    unsigned cell_containing_node_mi = cell_store[current_boundary_node_ci];
 
    for(unsigned current_boundary_node_local=1;current_boundary_node_local<no_boundary_nodes;current_boundary_node_local++){

        if(cells_per_node[current_boundary_node_ci]==2){
            
            for(unsigned i=0;i<2;i++){
                
                unsigned current_cell_containing_node_mi = (unsigned)cell_store[current_boundary_node_ci+i*array_sizes];
                
                if(current_cell_containing_node_mi!=cell_containing_node_mi){
                    cell_containing_node_mi = current_cell_containing_node_mi;
                    break;
                }
            }
        }
            
        unsigned cell_containing_node_ci = cell_containing_node_mi-1;
        
        mxArray* mx_cell_nodes = mxGetCell(prhs[1], cell_containing_node_ci);
        double* cell_nodes = mxGetPr(mx_cell_nodes);
        unsigned no_cell_nodes = mxGetN(mx_cell_nodes);
        
        for(unsigned current_cell_node_local=0;current_cell_node_local<no_cell_nodes;current_cell_node_local++){
            if(cell_nodes[current_cell_node_local]==current_boundary_node_mi){
                current_boundary_node_mi = cell_nodes[(current_cell_node_local+1)%no_cell_nodes];
                break;
            }
        }
        
        current_boundary_node_ci = current_boundary_node_mi-1;
        boundary_element[current_boundary_node_local] = current_boundary_node_mi;
    }
}

