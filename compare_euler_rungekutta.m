disp('busy');clear all; close all;

euler = load('euler');
rungekutta = load('rungekutta');

figure('outerposition',[100 100 700 700])
cellfun(@(x)patch(euler.vertices.position(x,1),euler.vertices.position(x,2),...
	[255,30,30]/255,'linewidth',2,'FaceAlpha',0,'edgecolor','g'),euler.cells.vertices);
hold on
cellfun(@(x)patch(rungekutta.vertices.position(x,1),rungekutta.vertices.position(x,2),...
	[255,30,30]/255,'linewidth',2,'FaceAlpha',0,'edgecolor','r'),rungekutta.cells.vertices);
axis equal

sum(sum((euler.vertices.position-rungekutta.vertices.position).^2,2))

