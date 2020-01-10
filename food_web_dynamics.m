function[output, solution] = food_web_dynamics(network_struct, metabolics_struct, parameter_set, sim_ID, display_figures)
% FOOD_WEB_DYNAMICS Use the set of differential equations specified in the
% Methods to simulate the dynamics of species abundances changing over
% time. This function is ONLY FOR N == S where N is the number of state 
% variables and S is the number of species, i.e. when there are NO plants 
% with pollinators or rewards. The FW dynamic model corresponds to running 
% traditional Allometric Trophic Network (ATN) dynamics, implemented
% following Martinez et al. 2012 Sustaining economic exploitation of 
% complex ecosystems in computational models of coupled human-natural 
% networks.
% It requires the following steps:
% (1) Set Run Identifiers: strings to indicate the network, metabolics, 
%       dynamics, and parameter set that were used
% (2) Format Inputs & Declare Parameters
% (3) Run the Simulation: uses the ode15s solver which various options for
%       events functions
% (4) Process Outputs: display a timeseries, save 'output' and/or
%       'solution'
%
%   output is a structure including .run_indentifiers, .final_biomasses, 
%       .persistence (species persistence), and .event_triggers like
%       species extinctions
%   solution is an N x timesteps array where each column corresponds to a 
%       timestep (starting at 0) and each row corresponds to a state 
%       variable (i.e. species or rewards) whose biomass is changing over 
%       time
%   network_struct is a structure including initial diversity (.S) &
%       integer indices of all guilds (e.g. .carnivores, .omnivores, etc.)
%   metabolics_struct is a structure including metabolic rate (.x, == 0 
%       for plants) & short-weighted trophic level (.swTL. == 1 for plants) 
%       for all species
%   parameter_set is a structure including all relevant (non-allometrically-
%       derived) parameter values
%   sim_ID is the (string) unique identifier for this simulation; when 
%       runnning in parallel, sim_ID is simply set to the integer simulation
%       run with padded 0s and all other details are saved in the 
%       run_identifiers (stored in the output)
%   display_figures is a logical flag, when == 1 plots a timeseries
%
% CITE THIS CODE AS FOLLOWS:
% Hale, K.R.S. (2020). Pollinators in food webs?Mutualistic interactions 
%   increase diversity, stability, and function in multiplex networks

% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%

% (1) Set Run Identifiers

file_ID = network_struct.fileID; % which niche model & pollination network used to construct network
mets_algorithm = metabolics_struct.mets_algorithm; % which algorithm used to create metabolic rates
web_strategy = metabolics_struct.web_strategy; % which consumer strategy used for metabolics (e.g. invertebrates, vertebrates)
parameter_ID = parameter_set.name; % which set of parameters

run_identifiers = struct('sim_ID',sim_ID,'file_ID',file_ID,'dynamics',...
    'food_web','mets_algorithm',mets_algorithm,'web_strategy',...
    web_strategy,'parameter_set',parameter_ID);

% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%

% (2) Format Inputs & Declare Parameters

% (2-1) Unpack interaction matrix and species roles from network_struct
% S: number of species
S = network_struct.S;

% I: interaction matrix, size is S x S
I = network_struct.I;

% species roles
plants = network_struct.wind; % here, all plants are wind-pollinated % pollinator-independent plants: wind-pollinated/selfing etc. 
consumers = setdiff(1:S,plants)';

% (2-2) Set allometric parameters
% x: allometrically-scaled mass-specific metabolic rate
x = metabolics_struct.x;

% r: allometrically-scaled maximum mass-specific growth rate of plants
r_wind = parameter_set.r_wind;

% yij: allometrically-scaled maximum consumsumption rate of i eating j
yij = 10*I; % AAAI 2012 invertebrates % 8*I; % Brose et al. 2006 invertebrates %

% eij: assimilation efficiency for eating organisms or floral resources
%   calculate the reciprocal directly since it is always used
eij_eating_vegetation = 1/0.66; % Martinez et al. 2012 % or % 1/0.45; % Brose et al. 2006 %
eij_eating_animals = 1/0.85; % Brose et al. 2006 & Martinez et al. 2012

% (2-4) Set other consumer-resource (ATN) parameters
% h: Hill coefficient (in functional response)
h = parameter_set.h; % 1 % Holling Type II % 2 % Holling Type III % 1.2 % Martinez et al. 2006 %

% B0: half-saturation density for i eating j (in functional response)
%   see discussion in SI and definition in Boit et al. 2012 
%   calculate B0^h directly since it is always used
B0 = parameter_set.B0;
B0 = I .* B0;
% B0(51:end,plants) = I(51:end,plants) * 30;
B0_h = B0.^h;

% omega: preference of i for eating j a.k.a. the generalist model 
%   see Williams 2008 for weak vs. strong generalist model
preference_model = parameter_set.preference_model;
omega = calc_preference(preference_model,I);

% competition_model: form of plant competition
competition_model = parameter_set.competition_model;

% K: plant-community-wide carrying capacity
K = parameter_set.K;

% (2-5) Set simulation parameters 
% initial_biomass: biomasses the simulation starts at
initial_biomass = parameter_set.initial_biomass;
initial_biomass = ones(S,1) * initial_biomass;

% num_timesteps: how long to run the simulation for
num_timesteps = parameter_set.time_steps;

% extinction_threshold: biomass at which species i is considered "extinct"
extinction_threshold = parameter_set.extinction_threshold; % (B(i) := 0, performed in event function, set options function for 'extinction')

% timer: length of time the simulation is allowed to run
timer = parameter_set.timer;  % set options function for 'stop_solve'

% (2-6) Allocate memory for simulation variables
B_h = initial_biomass.^h;
F = zeros(S,S);
consumption_rate = zeros(S,S);
extinct = zeros(S,1);

% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%

% (3) Run the Simulation

% (3-1) Define nested function for solving the differential equations
function [dBdt] = dynamics_RHS(~, B)
        
    B(B <= extinction_threshold) = 0;

    % dBdt is the rate of change for the biomass of species 1:S
    B_h = B.^h;

    % F is the matrix of the functional response for i eating j
    for j = 1:S
        for i = 1:S 
            if I(i,j) == 1
                F(i,j) = (omega(i,j) * B_h(j))/(((I(i,:) .* omega(i,:)) * B_h) + B0_h(i,j));
            end
        end
    end

    % plant competition
    competition = calc_competition(competition_model,B,plants,K,0); 
    
    % consumption_rate is the matrix of the rate at which i consumes j
    for j = 1:S
        for i = 1:S
            if I(i,j) == 1
                consumption_rate(i,j) = x(i) * yij(i,j) * B(i) * F(i,j);
            end
        end
    end
    
    % ASSEMBLING THE DIFFERENTIAL EQUATION VECTOR
    dBdt = zeros(S,1);

    for i = 1:numel(plants)
        dBdt(plants(i)) = competition * r_wind * B(plants(i))  - eij_eating_vegetation * sum(consumption_rate(:,plants(i)));
    end

    for i = 1:numel(consumers)
        dBdt(consumers(i)) = sum(consumption_rate(consumers(i),:)) - eij_eating_animals * sum(consumption_rate(:,consumers(i))) - x(consumers(i)) * B(consumers(i));
    end
    
    if any(isnan(dBdt))
        error('Rate of change is not a number!\n')
    end
    
    if ~isreal(dBdt)
        fprintf('Non-real derivative.\n')
    end
     
    dBdt(B <= extinction_threshold) = 0;
    
    end

% (3-2-1) Define an event function to stop the solver if simulation is 
% running too long. Make sure to start and stop the timer (3-5) & (3-7).
% Events stop the integration when value == 0 and isterminal == 1. 
function [value, isterminal, direction] = stop_solve(~,~)
    
    % Define the countdown timer in seconds
    % 	declared in simulation parameters (2-5)
    
    % The solver runs until this value is negative
    value = toc - timer;
    
    % The function should terminate solver execution
    isterminal = 1;
    
    % The direction does not matter
    direction = 0;
    
end

% (3-2-2) Define an event function to stop the solver when a species goes 
% extinct, allowing a hard set of its biomass to 0 (performed in (3-6)).
function [value, isterminal, direction] = extinction(~,B)
    
    % Define the threshold at which species are considered extinct
    % 	declared in simulation parameters (2-5)
    value = zeros(S,1);

    % The solver runs until this value is negative
    value(B <= extinction_threshold) = -1;
    
    % The function should terminate solver execution
    isterminal = ones(S,1)-extinct;
    
    % The direction does not matter
    direction = zeros(S,1);
    
end

% (3-3) Set the options including event function (can only use one event at
% a time)
options = odeset(...
    'RelTol',1e-7,... % Error tolerances should be modified with extinction threshold
    'AbsTol',1e-9,...
    'NonNegative',1:S,...  
    'Events',@extinction); 
    %'Events',@stop_solve);

% % (3-5) Start timer
% tic

% (3-6) Run the solver, continuing until num_timesteps has been reached (if
% extinctions events occur first)
solution_struct = ode45(@dynamics_RHS,[0 num_timesteps],initial_biomass,options);

event_time_step = solution_struct.xe;
biomasses_at_extinction = solution_struct.y(:,end);
extinct = biomasses_at_extinction < extinction_threshold;
biomasses_at_extinction(extinct) = 0;

while event_time_step < num_timesteps % && toc < timer
    solution_struct = odextend(solution_struct,[],num_timesteps,biomasses_at_extinction);
    event_time_step = solution_struct.x(end);
    
    biomasses_at_extinction = solution_struct.y(:,end);
    extinct = biomasses_at_extinction <= extinction_threshold;
    biomasses_at_extinction(extinct) = 0;
end

% (3-7) Save length of time the solver ran (in s)
% simulation_timed = toc;

% (3-8) Filter down solution struct to contain evaluations only at every
% time step
solution = deval(solution_struct,(0:num_timesteps));

% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%

% (4) Process Outputs

% (4-1) Save the final biomass distribution & last successful timestep 
final_biomasses = solution_struct.y(:,end);
final_timestep = solution_struct.x(:,end);
persistence = sum(final_biomasses(1:S) > extinction_threshold)/S;
 
% (4-2) Save the timesteps and biomasses at each extinction, conveniently
% recorded by the solver extinctions correspond to events. 
% NOTE: the timings of events and subsequent extinctions are of fuzzy
% accuracy for the following reasons:
%   (1) If more than one terminal event occurs during the first step, then 
%       only the first event registers and the solver continues integrating.
%   (2) Zeros are determined by sign crossings between steps. Therefore, 
%       zeros with an even number of crossings between steps can be missed.
% In other words, the species that triggered the extinction event might not
% actually go extinct right then, and thus may trigger another extinction
% event later. Also, two species may go extinct at the same timestep, and 
% the solver will only record the (first) one that triggered the event.
% This should not effect the hard setting of species' biomass to 0. I am
% thus recording the species that triggered extinction events for the sake
% of debugging, but extinction order should be deduced from the
% extinction_biomasses. 
extinction_biomasses = solution_struct.ye;
extinction_timesteps = solution_struct.xe;
event_triggers = solution_struct.ie;

% (4-3) Display timeseries figures if desired
if display_figures
    plot_timeseries(network_struct,solution);
end

% (4-5) Format  outputs for memory constraints on HPC cluster
output = struct('run_identifiers',run_identifiers,'S',S,'final_timestep',...
    final_timestep,'final_biomasses',final_biomasses,'extinction_timesteps',...
    extinction_timesteps,'extinction_biomasses',extinction_biomasses,...
    'persistence',persistence,'event_triggers',event_triggers);

% (4-6) Save outputs, only specifier is held in sim_ID
output_filename = strcat('output',sim_ID,'.mat');
save(output_filename, 'output');

solution_filename = strcat('solution',sim_ID,'.mat');
save(solution_filename, 'solution');
end
