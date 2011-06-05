disp('busy');clear all;close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_iterations = 500;
delta_t = 1e-3;

N = 10;
diffusion_speed = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

movie_logical = 1;

axis_values = [0 1 -1 1 0 1];
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

Dpp = cos(pi*node_positions(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

screen_size = get(0,'screensize');
figure_position =...
	[(screen_size(3)-5/4*screen_size(4))/2 0 5/4*screen_size(4) screen_size(4)];

hf = figure('outerposition',figure_position,'color','white');
% ha = axes('outerposition',[0.3 0.1 0.4 0.8]);

if ~movie_logical
    subplot(1,2,1)
end    

% trisurf(FEM_elements,node_positions(:,1),node_positions(:,2),Dpp);
% axis(axis_values);
% axis off;
% grid off;
% view(FEM_angle);
% drawnow;

x_true = linspace(0,1,1000);
y_true = cos(pi*x_true);

plot(x_true,y_true);
hold on;
plot(linspace(xy_min,xy_max,N),Dpp(node_positions(:,2)>(1-1/(2*N))),'r')
axis(axis_values);
xlabel('x');
ylabel('Concentration')
legend('Analytic','FEM')
title(['N = ',num2str(N)])
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iteration = 1:total_iterations

    previous_node_positions = node_positions;    

    [I,J,AV,MV] =...
        Stiff2D(delta_t,FEM_elements,node_positions,previous_node_positions);

    lap_matrix = sparse(I,J,AV);
    mass_matrix = sparse(I,J,MV);
        
    COEFF_MAT = mass_matrix + 0.5*delta_t*(diffusion_speed*lap_matrix);
    rhs = (mass_matrix - 0.5*delta_t*(diffusion_speed*lap_matrix))*Dpp;
        
    Dpp = COEFF_MAT\rhs;
    
    if movie_logical && ~rem(iteration,movie_period)

        cla;
        y_true = cos(pi*x_true)*exp(-pi^2*diffusion_speed*(delta_t*iteration));
        plot(x_true,y_true);
        hold on;
        plot(linspace(xy_min,xy_max,N),Dpp(node_positions(:,2)>(1-1/(2*N))),'r');
        axis(axis_values);
        drawnow;
    
    end
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~movie_logical
    
    subplot(1,2,2)
    y_true = cos(pi*x_true)*exp(-pi^2*(delta_t*iteration));
    plot(x_true,y_true);
    hold on;
    plot(linspace(xy_min,xy_max,N),Dpp(node_positions(:,2)>(1-1/(2*N))),'r');
    axis([0 1 0 1]);
    drawnow;

end