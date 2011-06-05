disp('busy');clear all;close all;

Dpp_view = [0 75];
temp_axis_values = 1.2*[-1 1 -1 1];
temp_axis_values_Dpp = [temp_axis_values 0 0.55];
caxis_vals = [0 0.6];
green1 = [50,180,50]/255;
green2 = [50,255,50]/255;
white = [255,255,0]/255;
colormap_val = [linspace(green1(1),white(1),300)' ...
    linspace(green1(2),white(2),300)' linspace(green1(3),white(3),300)'];
% colormap_val = 'summer';
% colormap_val = [0 150 256]/256;

alphabet = 'abcdefghijklmnopqrstuvwxyz';

folder_name = 'Dpp_gradient_lambda_2/';
saved_iterations = 8333:8333:24999;

if ~exist(['Figs/',folder_name],'dir')
	mkdir('Figs/',folder_name);
end

load(['Saves/',folder_name,'iteration_',num2str(saved_iterations(1)),'.mat'])
load(['Saves/',folder_name,'initial_save.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cell Image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fig_name = ['cells_',alphabet(1)];

% if exist(['Figs/',folder_name,fig_name,'.eps'],'file')
% 	yesno = input(['Are you sure you want to overwrite file ',fig_name,'.eps? '],'s');
% 	if ~numel(yesno) || ~(strcmp(yesno(1),'y'))
% 		error(['Could not delete file ',fig_name,'.eps'])
% 	end
% end

figure('outerposition',[100 100 600 600])%,'PaperPositionMode','auto')
axes('position',[0 0 1 1])
figure_loop(cells,node_positions,1,'r',1);
axis(temp_axis_values);
set(gca,'xcolor','w')
set(gca,'ycolor','w')
set(gca,'zcolor','w')
set(gca,'ticklength',[0 0])
axis equal

% saveas(gcf,['Figs/',folder_name,fig_name,'.eps'],'psc2')
% close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dpp Image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fig_name = ['Dpp_',alphabet(1)];

% if exist(['Figs/',folder_name,fig_name,'.eps'],'file')
% 	yesno = input(['Are you sure you want to overwrite file ',fig_name,'.eps? '],'s');
% 	if ~numel(yesno) || ~(strcmp(yesno(1),'y'))
% 		error(['Could not delete file ',fig_name,'.eps'])
% 	end
% end

figure('outerposition',[100 100 600 600])%,'PaperPositionMode','auto')
axes('position',[0 0 1 1])
trisurf(FEM_elements,FEM_node_positions(:,1),FEM_node_positions(:,2),Dpp,'linewidth',1)
grid off;
shading(shading_style)
axis(temp_axis_values_Dpp);
caxis(caxis_vals)
set(gca,'xcolor','w')
set(gca,'ycolor','w')
set(gca,'zcolor','w')
set(gca,'ticklength',[0 0])
view(Dpp_view);
colormap(colormap_val)

% saveas(gcf,['Figs/',folder_name,fig_name,'.eps'],'psc2')
% close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for unused_variable = 1:length(saved_iterations)
	
	if ~exist(['Saves/',folder_name,'iteration_',num2str(saved_iterations(unused_variable)),'.mat'],'file')
		error(['No save for iteration_',num2str(saved_iterations(unused_variable)),'.mat'])
	end
	
	load(['Saves/',folder_name,'iteration_',num2str(saved_iterations(unused_variable)),'.mat']);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
% 	fig_name = ['cells_',alphabet(unused_variable+1)];
	
% 	if exist(['Figs/',folder_name,fig_name,'.eps'],'file')
% 		yesno = input(['Are you sure you want to overwrite file ',fig_name,'.eps? '],'s');
% 		if ~numel(yesno) || ~(strcmp(yesno(1),'y'))
% 			error(['Could not delete file ',fig_name,'.eps'])
% 		end
% 	end
	
    figure('outerposition',[100 100 600 600])%,'PaperPositionMode','auto')
    axes('position',[0 0 1 1])
    figure_loop(cells,node_positions,1,'r',1);
    axis(temp_axis_values);
    set(gca,'xcolor','w')
    set(gca,'ycolor','w')
    set(gca,'zcolor','w')
    set(gca,'ticklength',[0 0])
    axis equal
	
% 	saveas(gcf,['Figs/',folder_name,fig_name,'.eps'],'psc2')
% 	close;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	fig_name = ['Dpp_',alphabet(unused_variable+1)];
% 	if exist(['Figs/',folder_name,fig_name,'.eps'],'file')
% 		yesno = input(['Are you sure you want to overwrite file ',fig_name,'.eps? '],'s');
% 		if ~numel(yesno) || ~(strcmp(yesno(1),'y'))
% 			error(['Could not delete file ',fig_name,'.eps'])
% 		end
% 	end
	
%     FEM_elements = FEM_elements(FEM_elements(:,1)>0);
	figure('outerposition',[100 100 600 600])%,'PaperPositionMode','auto')
	axes('position',[0 0 1 1])
	trisurf(FEM_elements,FEM_node_positions(:,1),FEM_node_positions(:,2),Dpp,'linewidth',1)
	grid off;
	shading(shading_style)
	axis(temp_axis_values_Dpp);
	caxis(caxis_vals)
    set(gca,'xcolor','w')
    set(gca,'ycolor','w')
    set(gca,'zcolor','w')
    set(gca,'ticklength',[0 0])
	view(Dpp_view);
	colormap(colormap_val)
	
% 	saveas(gcf,['Figs/',folder_name,fig_name,'.eps'],'psc2')
% 	close;
	
end