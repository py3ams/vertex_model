disp('busy');clear all; close all;

red1 = [255,30,30]/255;
green1 = [50,180,50]/255;
green2 = [50,255,50]/255;

figure('outerposition',[100 100 400 400],'PaperPositionMode','auto')
axes('position',[0.05 0.05 0.9 0.9])

cos60 = cos(pi/3);
sin60 = sin(pi/3);
% rt2 = 1/2;

vertex_positions = [0,sin60;cos60,2*sin60;1+cos60,2*sin60;2*cos60+1,sin60;1+cos60,0;cos60,0];%+0.2*(rand(6,2)-0.5);

angle_of_mitosis = pi/4.4;
cell_to_divide_centroid = mean(vertex_positions);

mitosis_line_gradient = tan(angle_of_mitosis);
mitosis_line_intercept =...
    cell_to_divide_centroid(2)-mitosis_line_gradient*cell_to_divide_centroid(1);

intersection_counter = 0;

no_cell_nodes = length(vertex_positions);

for current_node_local = 1:no_cell_nodes
    
    clockwise_node_local = mod(current_node_local,no_cell_nodes)+1;
    
    current_node_position = vertex_positions(current_node_local,:);
    clockwise_node_position = vertex_positions(clockwise_node_local,:);

    current_edge_gradient =...
        (clockwise_node_position(2)-current_node_position(2))/...
        (clockwise_node_position(1)-current_node_position(1));
    
    current_edge_intercept =...
        current_node_position(2)-current_edge_gradient*current_node_position(1);

    x_intersection =...
        (current_edge_intercept-mitosis_line_intercept)/...
        (mitosis_line_gradient-current_edge_gradient);
    
    x_max = max(current_node_position(1),clockwise_node_position(1));
    x_min = min(current_node_position(1),clockwise_node_position(1));
    
    if x_intersection > x_min && x_intersection < x_max
        
        intersection_counter = intersection_counter+1;
        
        if intersection_counter == 1
            
            new_node_1_position =...
                [x_intersection mitosis_line_gradient*x_intersection+mitosis_line_intercept];
            
        elseif intersection_counter == 2
            
            new_node_2_position =...
                [x_intersection mitosis_line_gradient*x_intersection+mitosis_line_intercept];
            
            break
            
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

original_centre = mean(vertex_positions);
FEM_vertex_positions = [vertex_positions;original_centre];

patchAS(vertex_positions,red1,5)
% axis equal
axis([-0.2 2.15 -0.2 2.15])
% axis tight
axis off

hold on

plot(cell_to_divide_centroid(1),cell_to_divide_centroid(2),'kx','linewidth',4,'markersize',20)

noise = 0.2;

% x_1 = 0.25+noise*(rand-0.5);
% y_1 = -sin60/cos60*x_1+sin60;
% 
% x_2 = noise*(rand-0.5)+1.75;
% y_2 = -sin60/cos60*x_2+(sin60+sin60/cos60*(2*cos60+1));

x_1 = new_node_2_position(1);
y_1 = new_node_2_position(2);

x_2 = new_node_1_position(1);
y_2 = new_node_1_position(2);

% plot(x_1,y_1,'bo','linewidth',3)
% plot(x_2,y_2,'bo','linewidth',3)

FEM_elements = [1 2 7; 2 3 7; 3 4 7; 4 5 7; 5 6 7; 6 1 7];

figure('outerposition',[100 100 400 400],'PaperPositionMode','auto')
axes('position',[0.05 0.05 0.9 0.9])

patch(vertex_positions(:,1),vertex_positions(:,2),red1,'linewidth',4,'FaceAlpha',0)
hold on

for i = 1:length(FEM_elements)
    patchAS(FEM_vertex_positions(FEM_elements(i,:),:),green1,1.5)
end
axis([-0.2 2.15 -0.2 2.15])
% axis equal
axis off

text(FEM_vertex_positions(1,1)-0.125,FEM_vertex_positions(1,2),'1','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(2,1)-0.025,FEM_vertex_positions(2,2)+0.1,'2','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(3,1)-0.025,FEM_vertex_positions(3,2)+0.1,'3','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(4,1)+0.05,FEM_vertex_positions(4,2),'4','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(5,1)-0.035,FEM_vertex_positions(5,2)-0.125,'5','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(6,1)-0.035,FEM_vertex_positions(6,2)-0.125,'6','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(7,1)-0.0375,FEM_vertex_positions(7,2)-0.175,'7','FontName',...
    'FixedWidth','FontSize',20)
    
vertex_positions = [vertex_positions;x_1 y_1;x_2 y_2];
cells = {[1 2 3 8 7],[8 4 5 6 7]};

figure('outerposition',[100 100 400 400],'PaperPositionMode','auto','color','white')
axes('position',[0.05 0.05 0.9 0.9])
patchAS(vertex_positions(cells{1},:),red1,5)
patchAS(vertex_positions(cells{2},:),red1,5)

axis([-0.2 2.15 -0.2 2.15])
axis off

FEM_vertex_positions = [vertex_positions;mean(vertex_positions(cells{1},:));mean(vertex_positions(cells{2},:))];
FEM_elements = [7 1 9; 1 2 9; 2 3 9; 3 8 9; 8 7 9; 8 4 10; 4 5 10; 5 6 10; 6 7 10; 7 8 10];

figure('outerposition',[100 100 400 400],'PaperPositionMode','auto','color','white')
axes('position',[0.05 0.05 0.9 0.9])

hold on
patch(vertex_positions(cells{1},1),vertex_positions(cells{1},2),red1,'linewidth',5,'FaceAlpha',0)
patch(vertex_positions(cells{2},1),vertex_positions(cells{2},2),red1,'linewidth',5,'FaceAlpha',0)

for i = 1:length(FEM_elements)
    patchAS(FEM_vertex_positions(FEM_elements(i,:),:),green1,1.5)
end
axis([-0.2 2.15 -0.2 2.15])
axis off

text(FEM_vertex_positions(1,1)-0.125,FEM_vertex_positions(1,2),'1','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(2,1)-0.025,FEM_vertex_positions(2,2)+0.1,'2','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(3,1),FEM_vertex_positions(3,2)+0.1,'3','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(4,1)+0.05,FEM_vertex_positions(4,2),'4','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(5,1),FEM_vertex_positions(5,2)-0.11,'5','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(6,1)-0.0225,FEM_vertex_positions(6,2)-0.125,'6','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(7,1)-0.1,FEM_vertex_positions(7,2)-0.075,'8','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(8,1)+0.05,FEM_vertex_positions(8,2)+0.05,'9','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(9,1),FEM_vertex_positions(9,2)-0.1,'10','FontName',...
    'FixedWidth','FontSize',20)
text(FEM_vertex_positions(10,1)+0.1,FEM_vertex_positions(10,2)-0.038,'11','FontName',...
    'FixedWidth','FontSize',20)

FEM_vertex_positions = [FEM_vertex_positions;original_centre];
FEM_elements1 = [7 1 11; 1 2 11; 2 3 11; 3 8 11;];
FEM_elements2 = [8 4 11; 4 5 11; 5 6 11; 6 7 11];

figure('outerposition',[100 100 400 400],'PaperPositionMode','auto')
axes('position',[0.05 0.05 0.9 0.9])
hold on
patch(vertex_positions(cells{1},1),vertex_positions(cells{1},2),red1,'linewidth',5,'FaceAlpha',0)
patch(vertex_positions(cells{2},1),vertex_positions(cells{2},2),red1,'linewidth',5,'FaceAlpha',0)
for i = 1:length(FEM_elements1)
    patchAS(FEM_vertex_positions(FEM_elements1(i,:),:),green1,1.5)
end
axis([-0.2 2.15 -0.2 2.15])
axis off

for i = 1:length(FEM_elements2)
    patchAS(FEM_vertex_positions(FEM_elements2(i,:),:),green2,1.5)
end
axis([-0.2 2.15 -0.2 2.15])
axis off


