function cell_figs_from_saves()

disp('busy');close all;

FEM_logical = 0;
refinement_level = 2;

save_figs_logical = 1;
initial_fig_logical = 1;

folder_name = 'ode_comparison/true_solution2/';
saved_iterations = [1600 16000];

% we set these outside function so both cell_fig and fem_fig have access to them
temp_axis_values = 0.58*[-1 1 -1 1];
% temp_axis_values = [-0.1 0.1 -0.07 0.13];
% apoptosis_figs
% temp_axis_values = [0.115 0.215 0.195 0.295];

if ~exist(['Figs/',folder_name],'dir')
    mkdir('Figs/',folder_name);
end

if initial_fig_logical
    
    % Cell fig
    
    load(['Saves/',folder_name,'initial_save.mat'])
    
    fig_name = ['Cells_0'];
    cell_vertices = cells.vertices;
    vertex_positions = vertices.position;
    cell_fig(cell_vertices,vertex_positions,fig_name,folder_name,save_figs_logical,temp_axis_values)
    
    % Dpp fig
    
    if FEM_logical
    
       % FEM_element_nodes = FEM_elements{refinement_level}.nodes;
       % FEM_node_positions = FEM_nodes{refinement_level}.position;
       % FEM_node_concentrations = FEM_nodes{refinement_level}.concentration;
       
       FEM_element_nodes = FEM_elements.nodes;
       FEM_node_positions = FEM_nodes.position;
       FEM_node_concentrations = FEM_nodes.concentration;
       
       fig_name = ['FEM_0'];
       fem_fig(FEM_element_nodes,FEM_node_positions,FEM_node_concentrations,...
          cell_vertices,vertex_positions,fig_name,folder_name,save_figs_logical,temp_axis_values)
       
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for unused_variable = 1:length(saved_iterations)
    
   % won't we get an error from the following line anyway?
    if ~exist(['Saves/',folder_name,'iteration_',num2str(saved_iterations(unused_variable)),'.mat'],'file')
        error(['No save for iteration_',num2str(saved_iterations(unused_variable)),'.mat'])
    end
    
    load(['Saves/',folder_name,'iteration_',num2str(saved_iterations(unused_variable)),'.mat']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cell_vertices = cells.vertices;
    vertex_positions = vertices.position;
    fig_name = ['Cells_',num2str(saved_iterations(unused_variable))];
    cell_fig(cell_vertices,vertex_positions,fig_name,folder_name,save_figs_logical,temp_axis_values);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    if FEM_logical
    
       %     FEM_element_nodes = FEM_elements{refinement_level}.nodes;
       %     FEM_node_positions = FEM_nodes{refinement_level}.position;
       %     FEM_node_concentrations = FEM_nodes{refinement_level}.concentration;
       
       FEM_element_nodes = FEM_elements.nodes;
       FEM_node_positions = FEM_nodes.position;
       FEM_node_concentrations = FEM_nodes.concentration;
       
       fig_name = ['FEM_',num2str(saved_iterations(unused_variable))];
       fem_fig(FEM_element_nodes,FEM_node_positions,FEM_node_concentrations,...
          cell_vertices,vertex_positions,fig_name,folder_name,save_figs_logical,temp_axis_values)
       
    end
    
end

end

% we delibrately pass fields rather than structures into the function to make the
% function more versatile
function cell_fig(cell_vertices,vertex_positions,fig_name,folder_name,save_figs_logical,temp_axis_values)

linewidth = 5;

figure('position',[100 100 500 500],'color','white','PaperPositionMode','auto')
axes('position',[0.01 0.01 0.98 0.98],'linewidth',2,'xcolor','white','ycolor','w',...
    'zcolor','w','ticklength',[0 0],'xtick',[],'ytick',[])

hold on
for current_cell = 1:length(cell_vertices)
    patchAS(vertex_positions(cell_vertices{current_cell},:),'r',linewidth)
end

axis(temp_axis_values)
% box on

if save_figs_logical
    saveas(gcf,['Figs/',folder_name,fig_name,'.eps'],'psc2')
end
    % close;

end

function fem_fig(FEM_element_nodes,FEM_node_positions,FEM_node_concentrations,...
    cell_vertices,vertex_positions,fig_name,folder_name,save_figs_logical,temp_axis_values)

cell_mesh_logical = true;
FEM_figs_logical = true;
Dpp_view = [0 90];
temp_axis_values_FEM = [temp_axis_values 0 0.1];
caxis_vals = [0 0.1];
green1 = [50,180,50]/255;
% green2 = [50,255,50]/255;
white = [255,255,0]/255;
colormap_val = [linspace(green1(1),white(1),300)' ...
    linspace(green1(2),white(2),300)' linspace(green1(3),white(3),300)'];
% colormap_val = 'summer';
% colormap_val = [0 150 256]/256;
shading_style = 'faceted';

if FEM_figs_logical
    
    figure('position',[100 100 500 500],'color','white','PaperPositionMode','auto')
    axes('position',[0.01 0.01 0.98 0.98],'linewidth',2,'xcolor','white','ycolor','w',...
        'zcolor','w','ticklength',[0 0],'xtick',[],'ytick',[])
    
    FEM_element_nodes = FEM_element_nodes(FEM_element_nodes(:,1)>0,:);
    
    if cell_mesh_logical

        trisurf(FEM_element_nodes,FEM_node_positions(:,1),FEM_node_positions(:,2),...
            zeros(length(FEM_node_positions(:,1)),1),'linewidth',1)
        
        grid off;axis off;shading(shading_style);axis(temp_axis_values);
        caxis(caxis_vals);colormap(colormap_val);view(Dpp_view);
        
        hold on
        cellfun(@(x)patch(vertex_positions(x,1),vertex_positions(x,2),'k','linewidth',3,'FaceAlpha',0),cell_vertices);
       
    else
        
        trisurf(FEM_element_nodes,FEM_node_positions(:,1),FEM_node_positions(:,2),...
            FEM_node_concentrations,'linewidth',1)
        
        grid off;axis off;shading(shading_style);axis(temp_axis_values_FEM);
        caxis(caxis_vals);colormap(colormap_val);view(Dpp_view);
        
    end
    
    if save_figs_logical
        saveas(gcf,['Figs/',folder_name,fig_name,'.eps'],'psc2')
    end
    % close;
    
end

end

