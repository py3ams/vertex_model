function cell_figs_from_saves()

disp('busy');close all;

cell_figs_logical = 0;
FEM_figs_logical = 1;

% 1- iteration number 2-letters
fig_labels = 2;

multiple_refinement_levels_logical = false;
refinement_level = 0;

save_figs_logical = 1;
initial_fig_logical = 1;

folder_name = 'drosophila_epidermis_unlimited_spi';
saved_iterations = [1000:1000:5000];

alphabet = 'abcdefghijklmnopqrstuvwxyz';
cell_alphabet_index = 1;
FEM_alphabet_index = 1;

% we set these outside function so both cell_fig and fem_fig have access to them
temp_axis_values = 1.0*[-1 1 -1 1];
% temp_axis_values = [-0.1 0.1 -0.07 0.13];
% apoptosis_figs
% temp_axis_values = [0.115 0.215 0.195 0.295];
% proliferation_figs
% temp_axis_values = [-0.125 0.025 -0.015 0.115];
% mesh figs (from realtime_refinement_comparison/true_iterations_8000_refinements_4)
% temp_axis_values = [-0.55 0.55 -0.52 0.58];
% T1 swaps
% temp_axis_values = [-0.3 -0.1 -0.4 -0.2];
% temp_axis_values = [-0.1 0.2 -0.02 0.22];

% temp_axis_values = [0.05 0.55 -0.25 0.25];

if save_figs_logical && ~exist(['Figs/',folder_name],'dir')
    mkdir('Figs/',folder_name);
end

if initial_fig_logical
    
    % Cell fig
    
    load(['Saves/',folder_name,'/initial_save.mat'])
    cell_vertices = cells.vertices;
    vertex_positions = vertices.position;
    cell_state = cells.state;
    
    if cell_figs_logical
        
       if fig_labels == 1
          fig_name = [folder_name,'_cells_0'];
       else
          fig_name = [folder_name,'_cells_',alphabet(cell_alphabet_index)];
          cell_alphabet_index = cell_alphabet_index+1;
       end
       
       cell_fig(cell_vertices,vertex_positions,fig_name,folder_name,save_figs_logical,temp_axis_values,cell_state)
        
    end
    
    % Dpp fig
    
    if FEM_figs_logical
    
       if multiple_refinement_levels_logical
       
          FEM_element_nodes = FEM_elements{refinement_level}.nodes;
          FEM_node_positions = FEM_nodes{refinement_level}.position;
          FEM_node_concentrations = FEM_nodes{refinement_level}.concentration;
          
       else
       
          FEM_element_nodes = FEM_elements.nodes;
          FEM_node_positions = FEM_nodes.position;
          FEM_node_concentrations = FEM_nodes.concentration;
          
       end
       
       if fig_labels == 1
          fig_name = [folder_name,'_FEM_0'];
       else
          fig_name = [folder_name,'_FEM_',alphabet(FEM_alphabet_index)];
          FEM_alphabet_index = FEM_alphabet_index+1;
       end
       
       fem_fig(FEM_element_nodes,FEM_node_positions,FEM_node_concentrations,...
          cell_vertices,vertex_positions,fig_name,folder_name,save_figs_logical,temp_axis_values)
       
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for unused_variable = 1:length(saved_iterations)
    
   % won't we get an error from the following line anyway?
    if ~exist(['Saves/',folder_name,'/iteration_',num2str(saved_iterations(unused_variable)),'.mat'],'file')
        error(['No save for iteration_',num2str(saved_iterations(unused_variable)),'.mat'])
    end
    
    load(['Saves/',folder_name,'/iteration_',num2str(saved_iterations(unused_variable)),'.mat']);
    cell_vertices = cells.vertices;
    vertex_positions = vertices.position;
    cell_state = cells.state;
    
    if cell_figs_logical
    
       if fig_labels == 1
          fig_name = [folder_name,'_cells_',num2str(saved_iterations(unused_variable))];
       else
          fig_name = [folder_name,'_cells_',alphabet(cell_alphabet_index)];
          cell_alphabet_index = cell_alphabet_index+1;
       end
   
        cell_fig(cell_vertices,vertex_positions,fig_name,folder_name,save_figs_logical,temp_axis_values,cell_state);
    
    end
        
    if FEM_figs_logical
    
       if multiple_refinement_levels_logical
       
          FEM_element_nodes = FEM_elements{refinement_level}.nodes;
          FEM_node_positions = FEM_nodes{refinement_level}.position;
          FEM_node_concentrations = FEM_nodes{refinement_level}.concentration;
       
       else
          
          FEM_element_nodes = FEM_elements.nodes;
          FEM_node_positions = FEM_nodes.position;
          FEM_node_concentrations = FEM_nodes.concentration;
          
       end
       
       if fig_labels == 1
          fig_name = [folder_name,'_FEM_',num2str(saved_iterations(unused_variable))];
       else
          fig_name = [folder_name,'_FEM_',alphabet(FEM_alphabet_index)];
          FEM_alphabet_index = FEM_alphabet_index+1;
       end
       fem_fig(FEM_element_nodes,FEM_node_positions,FEM_node_concentrations,...
          cell_vertices,vertex_positions,fig_name,folder_name,save_figs_logical,temp_axis_values)
       
    end
    
end

end

% we delibrately pass fields rather than structures into the function to make the
% function more versatile
function cell_fig(cell_vertices,vertex_positions,fig_name,folder_name,save_figs_logical,temp_axis_values,cell_state)

linewidth = 3;

figure('position',[100 100 500 500],'color','white','PaperPositionMode','auto')
axes('position',[0.01 0.01 0.98 0.98],'linewidth',2,'xcolor','white','ycolor','w',...
    'zcolor','w','ticklength',[0 0],'xtick',[],'ytick',[])

hold on
for current_cell = 1:length(cell_vertices)
   if cell_state(current_cell)==3
      patchAS(vertex_positions(cell_vertices{current_cell},:),'r',linewidth)
   else
      patchAS(vertex_positions(cell_vertices{current_cell},:),'r',linewidth)
   end
end

axis(temp_axis_values)

% xmin = 0.95;
% xmax = 1.45;
% ypos = -1.4;

xmin = 0.45;
xmax = 0.95;
ypos = -0.9;

plot([xmin xmax],[ypos ypos],'k','linewidth',4)
plot([xmin xmin],[ypos-0.05 ypos+0.05],'k','linewidth',4)
plot([xmax xmax],[ypos-0.05 ypos+0.05],'k','linewidth',4)
text(xmin+(xmax-xmin)/2,ypos+0.1,[num2str(xmax-xmin),' a.u.'],'FontSize',16,...
	'HorizontalAlignment','center','fontweight','bold','fontname','arial')

% box on

% if save_figs_logical
%    addpath('~/Documents/export_fig/')
%    export_fig(['Figs/',folder_name,'/',folder_name,'_',fig_name,'.eps'],'-nocrop');
% %    saveas(gcf,['Figs/',folder_name,fig_name,'.eps'],'psc2')
% end

if save_figs_logical
   addpath('~/Documents/export_fig/')
   export_fig(['Figs/',folder_name,'/',fig_name,'.eps'],'-nocrop');
end
    % close;

end

function fem_fig(FEM_element_nodes,FEM_node_positions,FEM_node_concentrations,...
    cell_vertices,vertex_positions,fig_name,folder_name,save_figs_logical,temp_axis_values)

cell_concentration_logical = true;

if cell_concentration_logical 
    Dpp_view = [0 60];
else
    Dpp_view = [0 90];
end
linewidth_cells = 5;
linewidth_elements = 1;
temp_axis_values_FEM = [temp_axis_values 0 50];
caxis_vals = [0 20];
green1 = [50,180,50]/255;
white1 = [255,255,0]/255;
colormap_val = [linspace(green1(1),white1(1),300)' ...
    linspace(green1(2),white1(2),300)' linspace(green1(3),white1(3),300)'];
shading_style = 'faceted';

figure('position',[100 100 500 500],'color','white','PaperPositionMode','auto','renderer','zbuffer')
axes('position',[0.01 0.01 0.98 0.98],'linewidth',2,'xcolor','white','ycolor','w',...
    'zcolor','w','ticklength',[0 0],'xtick',[],'ytick',[])

FEM_element_nodes = FEM_element_nodes(FEM_element_nodes(:,1)>0,:);

if ~cell_concentration_logical
    
    trisurf(FEM_element_nodes,FEM_node_positions(:,1),FEM_node_positions(:,2),...
        zeros(length(FEM_node_positions(:,1)),1),'linewidth',linewidth_elements)
    
    grid off;axis off;shading(shading_style);caxis(caxis_vals); view(Dpp_view);colormap(white)
    
    axis(temp_axis_values);
    
    cellfun(@(x)patch(vertex_positions(x,1),vertex_positions(x,2),'w','linewidth',linewidth_cells,'FaceAlpha',0),cell_vertices);
    
else
    
    trisurf(FEM_element_nodes,FEM_node_positions(:,1),FEM_node_positions(:,2),...
        FEM_node_concentrations,'linewidth',linewidth_elements)
    
    grid off;axis off;shading(shading_style);axis(temp_axis_values_FEM);
    caxis(caxis_vals);colormap(colormap_val);view(Dpp_view);
   
end

% if save_figs_logical
%    addpath('~/Documents/export_fig/')
%    export_fig(['Figs/',folder_name,'/',folder_name,'_',fig_name,'.eps'],'-nocrop');
% end
% % close;

if save_figs_logical
   addpath('~/Documents/export_fig/')
   export_fig(['Figs/',folder_name,'/',fig_name,'.eps'],'-nocrop');
end

end

