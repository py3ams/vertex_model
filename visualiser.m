function visualiser(cells,vertices,FEM_elements,FEM_nodes,axis_values,...
    axis_values_FEM,chemical_to_view,final_iteration_logical,iteration,...
    linewidth_cells,linewidth_elements,movie_logical,update_period,...
    view_FEM_concentration,view_FEM_mesh,view_initial_config,...
    view_iteration_number,view_number_cells)

global cell_figure_handle left_axes right_axes text_axes

if ~iteration
    
    if view_FEM_mesh
        
        if view_number_cells || view_iteration_number
            
            cell_figure_handle = figure('outerposition',[30 50 1200 715],'color','white');
            
            text_axes = axes('position',[0 0.9 1 0.05]);
            axis([0 1 0 1])
            axis off
            
            set_figure_text(cells,iteration,view_iteration_number,view_number_cells);
            
            left_axes = axes('position',[0 0 0.5 0.9]);
            right_axes = axes('position',[0.5 0 0.5 0.9]);
            
        else
            
            cell_figure_handle = figure('outerposition',[30 50 1200 650],'color','white');
            left_axes = axes('position',[0 0 0.5 1]);
            right_axes = axes('position',[0.5 0 0.5 1]);
            
        end
        
        set(cell_figure_handle,'CurrentAxes',left_axes);
        figure_loop(cells,vertices,linewidth_cells)
        axis(axis_values)
        axis off;
        
        set(cell_figure_handle,'CurrentAxes',right_axes);
        plot_FEM(cells,FEM_elements,FEM_nodes,axis_values_FEM,chemical_to_view,...
            linewidth_cells,linewidth_elements,vertices,view_FEM_concentration);
        
        
    else%if ~view_FEM_mesh
        
        if view_number_cells || view_iteration_number
            
            cell_figure_handle = figure('position',[300 50 650 715],'color','white');
            
            text_axes = axes('position',[0 0.9 1 0.05]);
            axis([0 1 0 1])
            axis off
            
            set_figure_text(cells,iteration,view_iteration_number,view_number_cells);
            
            left_axes = axes('position',[0 0 1 0.9]);
            
        else
            
            cell_figure_handle = figure('position',[300 50 650 650],'color','white');
            left_axes = axes('position',[0 0 1 1]);
            
        end
        
        figure_loop(cells,vertices,linewidth_cells)
        axis(axis_values)
        axis off;
        
    end
    
    if view_initial_config
        
        h1 = gcf;
        h2 = figure('position',get(gcf,'position'),'color','white');
        copyobj(get(h1,'children'),h2);
        
    end
    
elseif movie_logical && ~rem(iteration,update_period) || final_iteration_logical
    
    if view_number_cells || view_iteration_number
        
        set(cell_figure_handle,'CurrentAxes',text_axes);cla;
        set_figure_text(cells,iteration,view_iteration_number,view_number_cells)
        
    end
    
    set(cell_figure_handle,'CurrentAxes',left_axes);cla;
    figure_loop(cells,vertices,linewidth_cells)
    axis(axis_values)
    
    if view_FEM_mesh
        
        set(cell_figure_handle,'CurrentAxes',right_axes);cla;
        
        plot_FEM(cells,FEM_elements,FEM_nodes,axis_values_FEM,chemical_to_view,...
            linewidth_cells,linewidth_elements,vertices,view_FEM_concentration);
        
    end
    
    drawnow;
    
end

