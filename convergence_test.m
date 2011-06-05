disp('busy');tic;close all; clear all;

euler_truth  = load('convergence_test/delta_t_0001');
RK_truth = load('convergence_test/RK_delta_t_0001');

truth = euler_truth;

time_steps = [0.01;0.1;1;10;50;100];

approx = load('convergence_test/delta_t_001');
euler_approx(1) = approx.vertices;
% euler_difference(1) =...
% 	sqrt(sum(sum((approx.vertices.position - truth.vertices.position).^2)));
euler_difference(1) =...
	abs(approx.vertices.position(1) - truth.vertices.position(1));

approx = load('convergence_test/delta_t_01');
euler_approx(2) = approx.vertices;
% euler_difference(2) =...
% 	sqrt(sum(sum((approx.vertices.position - truth.vertices.position).^2)));
euler_difference(2) =...
	abs(approx.vertices.position(1) - truth.vertices.position(1));

approx = load('convergence_test/delta_t_1');
euler_approx(3) = approx.vertices;
% euler_difference(3) =...
% 	sqrt(sum(sum((approx.vertices.position - truth.vertices.position).^2)));
euler_difference(3) =...
	abs(approx.vertices.position(1) - truth.vertices.position(1));

approx = load('convergence_test/delta_t_10');
euler_approx(4) = approx.vertices;
% euler_difference(4) =...
% 	sqrt(sum(sum((approx.vertices.position - truth.vertices.position).^2)));
euler_difference(4) =...
	abs(approx.vertices.position(1) - truth.vertices.position(1));

approx = load('convergence_test/delta_t_50');
euler_approx(5) = approx.vertices;
% euler_difference(5) =...
% 	sqrt(sum(sum((approx.vertices.position - truth.vertices.position).^2)));
euler_difference(5) =...
	abs(approx.vertices.position(1) - truth.vertices.position(1));

approx = load('convergence_test/delta_t_100');
euler_approx(6) = approx.vertices;
% euler_difference(5) =...
% 	sqrt(sum(sum((approx.vertices.position - truth.vertices.position).^2)));
euler_difference(6) =...
	abs(approx.vertices.position(1) - truth.vertices.position(1));

figure('outerposition',[100 100 700 700])
semilogy(time_steps,euler_difference)
% plot(time_steps,euler_difference)
hold on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

truth = RK_truth;

approx = load('convergence_test/RK_delta_t_001');
RK_approx(1) = approx.vertices;
% RK_difference(1) =...
% 	sqrt(sum(sum((approx.vertices.position - truth.vertices.position).^2)));
RK_difference(1) =...
	abs(approx.vertices.position(1) - truth.vertices.position(1));

approx = load('convergence_test/RK_delta_t_01');
RK_approx(2) = approx.vertices;
% RK_difference(2) =...
% 	sqrt(sum(sum((approx.vertices.position - truth.vertices.position).^2)));
RK_difference(2) =...
	abs(approx.vertices.position(1) - truth.vertices.position(1));

approx = load('convergence_test/RK_delta_t_1');
RK_approx(3) = approx.vertices;
% RK_difference(3) =...
% 	sqrt(sum(sum((approx.vertices.position - truth.vertices.position).^2)));
RK_difference(3) =...
	abs(approx.vertices.position(1) - truth.vertices.position(1));

approx = load('convergence_test/RK_delta_t_10');
RK_approx(4) = approx.vertices;
% RK_difference(4) =...
% 	sqrt(sum(sum((approx.vertices.position - truth.vertices.position).^2)));
RK_difference(4) =...
	abs(approx.vertices.position(1) - truth.vertices.position(1));

approx = load('convergence_test/RK_delta_t_50');
RK_approx(5) = approx.vertices;
% RK_difference(4) =...
% 	sqrt(sum(sum((approx.vertices.position - truth.vertices.position).^2)));
RK_difference(5) =...
	abs(approx.vertices.position(1) - truth.vertices.position(1));

approx = load('convergence_test/RK_delta_t_100');
RK_approx(6) = approx.vertices;
% RK_difference(5) =...
% 	sqrt(sum(sum((approx.vertices.position - truth.vertices.position).^2)));
RK_difference(6) =...
	abs(approx.vertices.position(1) - truth.vertices.position(1));

semilogy(time_steps,RK_difference)
% plot(time_steps,RK_difference)
