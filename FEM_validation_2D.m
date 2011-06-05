disp('busy');clear all;close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_iterations = 100;
delta_t = 1e-3;

N = 20;
diffusion_speed = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axis_values = [0 1 0 1 -1 1];

FEM_angle = [0,30];
movie_period = 1;
frame_counter = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xy_min = 0;
xy_max = 1;

x = linspace(xy_min,xy_max,N);
y = linspace(xy_min,xy_max,N);

[x,y] = meshgrid(x,y);

x = x(:); y = y(:);

FEM_elements = delaunay(x,y);

node_positions = [x(:) y(:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dpp = cos(pi*node_positions(:,1)).*cos(pi*node_positions(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

screen_size = get(0,'screensize');
figure_position =...
	[(screen_size(3)-5/4*screen_size(4))/2 0 5/4*screen_size(4) screen_size(4)];

hf = figure('outerposition',figure_position,'color','white');
% ha = axes('outerposition',[0.3 0.1 0.4 0.8]);

subplot(1,2,1)
trisurf(FEM_elements,node_positions(:,1),node_positions(:,2),Dpp);
axis(axis_values);
axis off;
grid off;
view(FEM_angle);

subplot(1,2,2)
x_true = linspace(0,1,N);
y_true = linspace(0,1,N);

[x_true,y_true] = meshgrid(x_true,y_true);

z_true = cos(pi*x_true).*cos(pi*y_true);

surf(x_true,y_true,z_true);
axis(axis_values);
axis off;
grid off;
view(FEM_angle);
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iteration = 1:total_iterations

    previous_node_positions = node_positions;    

    [I,J,AV,MV] =...
        Stiff2D(delta_t,FEM_elements,node_positions,previous_node_positions);

    A = sparse(I,J,AV);
    M = sparse(I,J,MV);

    mass_matrix = M;
    lap_matrix = A;
        
    COEFF_MAT = mass_matrix + 0.5*delta_t*diffusion_speed*lap_matrix;
    rhs = (mass_matrix - 0.5*delta_t*diffusion_speed*lap_matrix)*Dpp;

%     COEFF_MAT = mass_matrix + delta_t*diffusion_speed*lap_matrix;
%     rhs = mass_matrix*Dpp;

    Dpp = COEFF_MAT\rhs;
    
    if ~rem(iteration,movie_period)
    
        subplot(1,2,1)
        cla;
        trisurf(FEM_elements,node_positions(:,1),node_positions(:,2),Dpp);
        axis(axis_values);
        axis off;
        grid off;
        view(FEM_angle);
        
        subplot(1,2,2)
        cla;
        z_true = cos(pi*x_true).*cos(pi*y_true).*exp(-2*pi^2*(delta_t*iteration));

        surf(x_true,y_true,z_true);
        axis(axis_values);
        axis off;
        grid off;
        view(FEM_angle);
        drawnow;
        
        frame_counter = frame_counter+1;
    
    end
        
    max_edge_length = sqrt(2*(mean(max(node_positions)-min(node_positions))/N).^2);
    max_node_speed =...
        max(sqrt(sum((node_positions-previous_node_positions).^2,2)))/delta_t;

end

N
sum(abs(Dpp-z_true(:)))/N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%