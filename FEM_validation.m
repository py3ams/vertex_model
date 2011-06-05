disp('busy');clear all;close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mex Stiff2D.cpp

total_iterations = 200;
delta_t = 1;

N = 40;

diffusion_speed = 0.0001;
% diffusion_speed = 1;
% source_magnitude = 0.1;
source_magnitude = 0;
source_radius = 0.2;
% degradation_constant = 0.01;
degradation_constant = 0;

mesh_growth_speed = 5e-4;
% mesh_growth_speed = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

movie_logical = 1;

% axis_values = [-0.5 0.5 -0.5 0.5 0 1];
axis_values = [-1.0 1.0 -1.0 1.0 0 2];
% axis_values = [0 1 0 1 -1 1];
% axis_values = 'equal';
FEM_angle = [0,45];
movie_period = 10;
frame_counter = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xy_min = -0.5;
xy_max = 0.5;

% xy_min = 0;
% xy_max = 1;

x = linspace(xy_min,xy_max,N);
y = linspace(xy_min,xy_max,N);

[x,y] = meshgrid(x,y);

x = x(:); y = y(:);

FEM_elements = delaunay(x,y);

node_positions = [x(:) y(:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dpp = zeros(N^2,1);

Dpp(sqrt(sum(node_positions.^2,2))<source_radius) = 1;
% Dpp = cos(pi*node_positions(:,1)).*cos(pi*node_positions(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_Dpp = zeros(total_iterations+1,1);
max_peclet_number = zeros(total_iterations,1);

total_Dpp(1) = CalculateTotalDpp(Dpp,FEM_elements,node_positions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

screen_size = get(0,'screensize');
figure_position =...
	[(screen_size(3)-5/4*screen_size(4))/2 0 5/4*screen_size(4) screen_size(4)];

hf = figure('outerposition',figure_position,'color','white');
% ha = axes('outerposition',[0.3 0.1 0.4 0.8]);

if ~movie_logical
    subplot(1,2,1)
end    

trisurf(FEM_elements,node_positions(:,1),node_positions(:,2),Dpp);
axis(axis_values);
axis off;
grid off;
view(FEM_angle);
drawnow;

% x_true = linspace(0,1,1000);
% y_true = cos(pi*x_true);
% 
% plot(x_true,y_true);
% hold on;
% plot(linspace(xy_min,xy_max,N),Dpp(node_positions(:,2)>(1-1/(2*N))),'r')
% axis([0 1 0 1]);
% drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% previous_node_positions = node_positions;    

for iteration = 1:total_iterations
        
    previous_node_positions = node_positions;    
    
    if iteration < total_iterations/2
        xy_min = xy_min-mesh_growth_speed;
        xy_max = xy_max+mesh_growth_speed;
    else
        xy_min = xy_min+mesh_growth_speed;
        xy_max = xy_max-mesh_growth_speed;
    end
    
    x = linspace(xy_min,xy_max,N);
    y = linspace(xy_min,xy_max,N);

    [x,y] = meshgrid(x,y);
    
    node_positions = [x(:) y(:)];

    [row_index,column_index,AV,MV,WV] =...
        Stiff2D(delta_t,FEM_elements,node_positions,previous_node_positions);

    A = sparse(row_index,column_index,AV);
    M = sparse(row_index,column_index,MV);
    W = sparse(row_index,column_index,WV);
    
    source_functions = zeros(length(node_positions),1);
    source_functions(sqrt(sum(node_positions.^2,2))<source_radius) =...
        source_magnitude;
               
%     COEFF_MAT = M + delta_t*(diffusion_speed*A);
    COEFF_MAT = M + delta_t*(diffusion_speed*A+W);
    
%     max(max(diffusion_speed*A))
%     max(max(W1))
    
    [row_index_prev,column_index_prev,MV_prev] = Stiff2DMonly(FEM_elements,previous_node_positions);
    M_prev = sparse(row_index_prev,column_index_prev,MV_prev);
    
    rhs = M_prev*Dpp*(1-degradation_constant) + delta_t*M*source_functions;
    
    Dpp = COEFF_MAT\rhs;
        
    if movie_logical && ~rem(iteration,movie_period)
    
        current_axis_values = axis;
        trisurf(FEM_elements,node_positions(:,1),node_positions(:,2),Dpp);
        axis(axis_values);
        axis off;
        grid off;
        view(FEM_angle);
        drawnow;
    
%         frame_counter = frame_counter+1;
%         F(frame_counter) = getframe(hf);

%         cla;
%         y_true = cos(pi*x_true)*exp(-pi^2*(delta_t*iteration));
%         plot(x_true,y_true);
%         hold on;
%         plot(linspace(xy_min,xy_max,N),Dpp(node_positions(:,2)>(1-1/(2*N))),'r');
%         axis([0 1 0 1]);
%         drawnow;
    
    end
        
    total_Dpp(iteration+1) =...
        CalculateTotalDpp(Dpp,FEM_elements,node_positions);
    
    max_edge_length = sqrt(2*(mean(max(node_positions)-min(node_positions))/N).^2);
    max_node_speed =...
        max(sqrt(sum((node_positions-previous_node_positions).^2,2)))/delta_t;
    
    max_peclet_number(iteration) =...
        max_edge_length*max_node_speed/diffusion_speed;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~movie_logical
    subplot(1,2,2)
    trisurf(FEM_elements,node_positions(:,1),node_positions(:,2),Dpp);
    axis(axis_values);
    axis off
    grid off;
    view(FEM_angle);
    drawnow;  
end

figure('outerposition',figure_position);

total_time = total_iterations*delta_t;
time_range = linspace(0,total_time,total_iterations+1);

subplot(1,2,1)
plot(time_range,total_Dpp)
xlabel('Time')
ylabel('Total Dpp')
axis([0 total_time 0.99*min(total_Dpp) 1.01*max(total_Dpp)])

subplot(1,2,2)
plot(time_range(2:end),max_peclet_number)
hold on
xlabel('Time')
ylabel('Maximum peclet number')

% for i=1:10
%     F(frame_counter+1:frame_counter+10) = getframe(gcf);
% end