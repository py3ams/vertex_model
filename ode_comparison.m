disp('busy');clear all;close all;

iterations_in_test_solution = 4000;
ode_solver_type = 2;

% simulation_name = 'ode_comparison/true_solution';
simulation_name = ['iterations_',num2str(iterations_in_test_solution),...
   '_ode_solver_',num2str(ode_solver_type)];

true_solution_initial_vars = load('Saves/ode_comparison/true_solution/initial_save');
test_solution_initial_vars = load(['Saves/',simulation_name,'/initial_save']);

assert(sum(sum(true_solution_initial_vars.vertices.position-test_solution_initial_vars.vertices.position))==0)

true_solution_vars = load('Saves/ode_comparison/true_solution/final_save');
test_solution_vars = load(['Saves/',simulation_name,'/final_save']);

solution_error = sqrt(sum(sum((true_solution_vars.vertices.position-test_solution_vars.vertices.position).^2,2)));

disp(['solution error : ',num2str(solution_error)])