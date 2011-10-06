disp('busy');close all; clear all;

load Saves/cellular_production_and_ingestion/final_save

statistics_counter = stats.counter;
total_concentration = stats.total_concentration(1:statistics_counter,:);
total_chemical_released = stats.chemical_source(1:statistics_counter,1);
total_internal_chemical = stats.internal_chemical_quantity(1:statistics_counter,5);
total_dead_chemical = stats.dead_chemical(1:statistics_counter);
net_chemical = total_chemical_released-total_internal_chemical-total_dead_chemical;

error_function = sqrt(1/total_time*delta_t*sum((1-net_chemical./total_concentration).^2));

disp(['Error : ',num2str(error_function)])
