% This script runs Rewards Only Food Web (FW) simulations. 
% Find & replace "RO" with "RP" for the Rewards Plus FW treatment.
% I recommend running this script in a folder to catch all 24,276
% output files and associated solution files. See "food_web_dynamics.m" 
% for detailed explanation the dynamics. 
%
% This code requires the parallel computing toolbox. To run in serial,
% follow the instructions in the script. With six cores on a MacBook with 
% 2.2 GHz Intel i7 Processor and 16 GB Memory, this code takes 
% approximately 3 hours to complete. 
% 
% CITE THIS CODE AS FOLLOWS:
% Hale, K.R.S. (2020). Pollinators in food webs: Mutualistic interactions 
%   increase diversity, stability, and function in multiplex networks

disp('Loading network structures and parameters...')

% Rewards Only (RO) Food Webs (FW)
load('networks_RO_FW.mat')
load('metabolics_RO_FW.mat')
load('simulation_parameters.mat','parameter_set')

networks = networks_RO_FW;
metabolics = metabolics_RO_FW;

parameter_set = parameter_set.food_web;

num_networks = numel(networks);
num_parameter_sets = numel(parameter_set);
num_simulations = num_networks * num_parameter_sets;

% don't display timeseries
display_figures = 0;

% comment out the following two lines to run in serial: 
numcores = feature('numcores');
parpool('local',numcores);

% use for instead of parfor to run in serial:
parfor i = 1:num_simulations
    this_network = mod(i,num_networks);
    this_sample = floor(i/(num_networks)) + 1;
    
    disp(['Running network: ', sprintf('%05d',i)])
    
    if (this_network == 0)
        this_network = num_networks;
        this_sample = this_sample - 1;
    end
        
    food_web_dynamics(networks(this_network), metabolics(this_network), ...
        parameter_set(this_sample), sprintf('%05d',i), display_figures);

end