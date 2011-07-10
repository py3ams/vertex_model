disp('busy');clear all;close all;

FEM_saves_logical = true;

Dpp_view = [0 75];
temp_axis_values = 0.6*[-1 1 -1 1];
temp_axis_values_FEM = [temp_axis_values 0 0.1];
caxis_vals = [0 0.1];
green1 = [50,180,50]/255;
green2 = [50,255,50]/255;
white = [255,255,0]/255;
colormap_val = [linspace(green1(1),white(1),300)' ...
    linspace(green1(2),white(2),300)' linspace(green1(3),white(3),300)'];
% colormap_val = 'summer';
% colormap_val = [0 150 256]/256;
shading_style = 'faceted';

refinement_level = 2;

folder_name = 'realtime_refinement_comparison/true_iterations_8000_refinements_4/';
saved_iterations = 800:800:8000;

if ~exist(['Figs/',folder_name],'dir')
    mkdir('Figs/',folder_name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['Saves/',folder_name,'initial_save.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cell Image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[200 200 500 500],'color','white','PaperPositionMode','auto')
axes('position',[0.01 0.01 0.98 0.98],'linewidth',2,'xcolor','white','ycolor','w',...
    'zcolor','w','ticklength',[0 0],'xtick',[],'ytick',[])

for current_cell = 1:length(cells.vertices)
    hold on
    patchAS(vertices.position(cells.vertices{current_cell},:),'r',2)
end

axis(temp_axis_values)
% box on

fig_name = ['Cells_initial'];
saveas(gcf,['Figs/',folder_name,fig_name,'.eps'],'psc2')
% close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dpp Image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FEM_saves_logical

    figure('position',[200 200 500 500],'color','white','PaperPositionMode','auto')
    axes('position',[0.01 0.01 0.98 0.98],'linewidth',2,'xcolor','white','ycolor','w',...
        'zcolor','w','ticklength',[0 0],'xtick',[],'ytick',[])
    
    trisurf(FEM_elements{refinement_level}.nodes(FEM_elements{refinement_level}.nodes(...
        :,1)>0,:),FEM_nodes{refinement_level}.position(:,1),FEM_nodes{refinement_level...
        }.position(:,2),FEM_nodes{refinement_level}.concentration,'linewidth',1)
    
    grid off;
    axis off;
    shading(shading_style);
    axis(temp_axis_values_FEM);
    caxis(caxis_vals);
    colormap(colormap_val);
    view(Dpp_view);
    
    fig_name = ['FEM_initial'];
    saveas(gcf,['Figs/',folder_name,fig_name,'.eps'],'psc2')
    % close;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for unused_variable = 1:length(saved_iterations)
    
    if ~exist(['Saves/',folder_name,'iteration_',num2str(saved_iterations(unused_variable)),'.mat'],'file')
        error(['No save for iteration_',num2str(saved_iterations(unused_variable)),'.mat'])
    end
    
    load(['Saves/',folder_name,'iteration_',num2str(saved_iterations(unused_variable)),'.mat']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure('position',[200 200 500 500],'color','white','PaperPositionMode','auto')
    axes('position',[0.01 0.01 0.98 0.98],'linewidth',2,'xcolor','white','ycolor','w',...
        'zcolor','w','ticklength',[0 0],'xtick',[],'ytick',[])
    
    for current_cell = 1:length(cells.vertices)
        hold on
        patchAS(vertices.position(cells.vertices{current_cell},:),'r',2)
    end
    
    %     box on
    
    axis(temp_axis_values)
    
    fig_name = ['Cells_',num2str(saved_iterations(unused_variable))];
    saveas(gcf,['Figs/',folder_name,fig_name,'.eps'],'psc2')
    %     close;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if FEM_saves_logical
        
        figure('position',[200 200 500 500],'color','white','PaperPositionMode','auto')
        axes('position',[0.01 0.01 0.98 0.98],'linewidth',2,'xcolor','white','ycolor','w',...
            'zcolor','w','ticklength',[0 0],'xtick',[],'ytick',[])
        
        trisurf(FEM_elements{refinement_level}.nodes(FEM_elements{refinement_level}.nodes(...
            :,1)>0,:),FEM_nodes{refinement_level}.position(:,1),FEM_nodes{refinement_level...
            }.position(:,2),FEM_nodes{refinement_level}.concentration,'linewidth',1)
        
        grid off;
        axis off;
        shading(shading_style);
        axis(temp_axis_values_FEM);
        caxis(caxis_vals);
        colormap(colormap_val);
        view(Dpp_view);
        
        fig_name = ['FEM_',num2str(saved_iterations(unused_variable))];
        saveas(gcf,['Figs/',folder_name,fig_name,'.eps'],'psc2')
        %     close;
    
    end
    
end