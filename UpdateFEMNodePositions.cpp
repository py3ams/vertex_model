#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if (nrhs != 4) {
        mexErrMsgTxt("4 input arguments required");
    }
    
    unsigned no_cells = mxGetM(prhs[0]);
    unsigned no_vertices = mxGetM(prhs[1]);
    unsigned no_FEM_nodes = mxGetM(prhs[3]);
    
    double* node_positions = mxGetPr(prhs[1]);
    double* cell_areas = mxGetPr(prhs[2]);
    double* edges = mxGetPr(prhs[3]);
    
    plhs[0] = mxCreateDoubleMatrix(no_FEM_nodes, 2, mxREAL);
    double* FEM_node_positions = mxGetPr(plhs[0]);
    
    // for the normal FEM node positions (i.e. not the edge nodes or the centroid nodes,
    // we can simply copy the values from previous position. if the node does not actually 
    // exist its position should be 0 so this will be copied across
    for(unsigned current_node_local=0;current_node_local<no_vertices;current_node_local++){
        
        for(unsigned dim=0;dim<2;dim++){
            
            FEM_node_positions[current_node_local+dim*no_FEM_nodes] =
                    node_positions[current_node_local+dim*no_vertices];
            
        }
    }
    
    unsigned FEM_index_zero_cell = no_FEM_nodes-no_cells;
    
    // for the edge nodes we must first check if the node actually exists, then work out its position
    // as the mid-point of the edge it is on.
    for(unsigned current_FEM_node_local=no_vertices;current_FEM_node_local<FEM_index_zero_cell;current_FEM_node_local++){
        
        unsigned edge_vertex_1_mi = (unsigned)edges[current_FEM_node_local];
        unsigned edge_vertex_2_mi = (unsigned)edges[current_FEM_node_local+no_FEM_nodes];
        
        if(edge_vertex_1_mi>0){
            
            unsigned edge_vertex_1_ci = edge_vertex_1_mi-1;
            unsigned edge_vertex_2_ci = edge_vertex_2_mi-1;
            
            FEM_node_positions[current_FEM_node_local] =
                    0.5*(node_positions[edge_vertex_1_ci]+
                    node_positions[edge_vertex_2_ci]);
            
            FEM_node_positions[current_FEM_node_local+no_FEM_nodes] =
                    0.5*(node_positions[edge_vertex_1_ci+no_vertices]+
                    node_positions[edge_vertex_2_ci+no_vertices]);
            
        }
    }
    
    // for the centroid nodes we work out the centroid of their cell. it is important not 
    // to attempt this if the cell does not exist, as it may create NaN's rather than 0's 
    // in the position.
    for(unsigned current_cell=0;current_cell<no_cells;current_cell++){
        
        mxArray* mx_cell_nodes = mxGetCell(prhs[0], current_cell);
        unsigned no_cell_nodes = mxGetN(mx_cell_nodes);
        
        if(no_cell_nodes>0){
            
            double* cell_nodes = mxGetPr(mx_cell_nodes);
            
            double cell_area = cell_areas[current_cell];
            
            double cx = 0, cy = 0;
            double x_i, x_iplus1, y_i, y_iplus1;
            
            for(unsigned i=0;i<no_cell_nodes;i++){
                
                unsigned current_cell_node_ci = (unsigned)cell_nodes[i]-1;
                unsigned clockwise_cell_node_ci = (unsigned)cell_nodes[(i+1)%no_cell_nodes]-1;
                
                x_i = node_positions[current_cell_node_ci];
                x_iplus1 = node_positions[clockwise_cell_node_ci];
                
                y_i = node_positions[current_cell_node_ci+no_vertices];
                y_iplus1 = node_positions[clockwise_cell_node_ci+no_vertices];
                
                cx = cx + (x_i + x_iplus1)*(x_i*y_iplus1 - x_iplus1*y_i);
                cy = cy + (y_i + y_iplus1)*(x_i*y_iplus1 - x_iplus1*y_i);
            }
            
            cx = cx/(6.0*(cell_area));
            cy = cy/(6.0*(cell_area));
            
            FEM_node_positions[FEM_index_zero_cell+current_cell] = -cx;
            FEM_node_positions[FEM_index_zero_cell+current_cell+no_FEM_nodes] = -cy;
            
        }
    }
}
