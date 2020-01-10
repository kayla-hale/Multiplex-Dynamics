function[production, consumption_on, consumption_by, metabolic_loss, realized_benefit] = calc_multiplex_flows(network_struct, metabolics_struct, parameter_set, B)
%CALC_MULTIPLEX_FLOWS This script is identical to "multiplex_dynamics.m"
%but is hard-wired to output flows (rates of change) at a specific biomass
%distribution B instead of integrating over timesteps to yield a 
%timeseries. As such, it requires only formatting inputs & parameters and
%calculating flows.
%
%   production, consumption_on, consumption_by, & metabolic_loss are double
%       vectors of those flows for each state variable. See Methods for a 
%       detailed definition of flows for each species & guild. 
%   realized_benefit is the outupt vector from the calc_pollination_benefit
%       function, corresponding to realized benefit from pollination to
%       each plant w/ pollinator from all its pollinators
%   network_struct is a structure including initial diversity (.S) &
%       integer indices of all guilds (e.g. .carnivores, .omnivores, etc.)
%   metabolics_struct is a structure including metabolic rate (.x, == 0 
%       for plants & rewards) & short-weighted trophic level (.swTL. == 1
%       for plants) for all species
%   parameter_set is a structure including all relevant (non-allometrically-
%       derived) parameter values
%   B is an Nx1 vector of biomasses, at which to evaluate flows, where N is
%       the number of state variables
%
% CITE THIS CODE AS FOLLOWS:
% Hale, K.R.S. (2019). Pollinators in food webs: Mutualistic interactions 
%   increase diversity, stability, and function in multiplex networks

% (2) Format Inputs & Declare Parameters

% (2-1) Get interaction matrix and species roles from network_struct
% S: number of species
S = network_struct.S;

% N: number of nodes (plants w/ pollinators have both vegetative (B) and
%   rewards (R) nodes)
N = network_struct.N;

% I: interaction matrix, size is N x N
I = network_struct.I;

% species roles (vectors of indices)
plants = network_struct.plants;
wind = network_struct.wind; % plants w/o pollinators (e.g. wind-pollinated, selfing etc.)
app = network_struct.app; % plants w/ pollinators' vegetative nodes (animal-pollinated plants [app])
rewards = network_struct.rewards; % plants w/ pollinators' floral rewards nodes
consumers = setdiff(1:S,plants)'; % all consumer guilds (e.g. herbivores, omnivores, etc.)
pollinators = network_struct.pollinators; % just pollinators

% (2-2) Set pollination parameters
% beta: production rate for floral rewards
beta = parameter_set.beta;

% s: self-limitation rate for floral rewards
s = parameter_set.s;

% kappa: cost of producing floral rewards
kappa = parameter_set.kappa;

% benefit_model: the form of the functional response of benefit accrual 
%   due to pollinator consumption of floral rewards
benefit_model = parameter_set.benefit_model;

% benefit_coefficient: shape parameter for functional response of benefit
%   accrual due to pollinator consumption of floral rewards
benefit_coefficient = parameter_set.benefit_coefficient;

% (2-3) Set allometric parameters
% x: allometrically-scaled mass-specific metabolic rate
x = metabolics_struct.x;

% r: allometrically-scaled maximum mass-specific growth rate of plants
r_wind = parameter_set.r_wind;
r_app = parameter_set.r_app;

% yij: allometrically-scaled maximum consumsumption rate of i eating j
yij = 10*I; % AAAI 2012 invertebrates % 8*IR; % Brose et al. 2006 invertebrates %

% eij: assimilation efficiency for eating organisms or floral resources
%   calculate the reciprocal directly since it is always used
eij_eating_vegetation = 1/0.66; % Martinez et al. 2012 % or % 1/0.45; % Brose et al. 2006 %
eij_eating_animals = 1/0.85; % Brose et al. 2006 & Martinez et al. 2012
eij_eating_rewards = 1; % just a guess

% (2-4) Set other consumer-resource (ATN) parameters
% h: Hill coefficient (in functional response)
h = parameter_set.h; % 1 % Holling Type II % 2 % Holling Type III % 1.2 % Martinez et al. 2006

% B0: half-saturation density for i eating j (in functional response)
%   see discussion in SI and definition in Boit et al. 2012 
%   calculate B0^h directly since it is always used
B0 = parameter_set.B0;
B0_floral = parameter_set.B0_floral;
B0 = ones(N,1) * B0;
B0(rewards) = B0_floral;
B0_h = B0.^h;

% omega: preference of i for eating j a.k.a. the generalist model 
%   see Williams 2008 for weak vs. strong generalist model
preference_model = parameter_set.preference_model;
preference_model_pollinators = parameter_set.preference_model_pollinators;
omega = calc_preference(preference_model,I);
omega(pollinators,:) = calc_preference(preference_model_pollinators,I(pollinators,:));

% competition_model: form of plant competition
competition_model = parameter_set.competition_model;

% K: plant-community-wide carrying capacity
K = parameter_set.K;

% extinction_threshold: biomass at which species i is considered "extinct"  
extinction_threshold = parameter_set.extinction_threshold; % (B(i) := 0, performed in event function, set options function for 'extinction')

% -----------------------------------------------------------------------%
% -----------------------------------------------------------------------%

% (3) Calculate the flows

B(B <= extinction_threshold) = 0;

% dBdt is the rate of change for the biomass of species 1:S and app's
% floral resources, indexed as S+1:N
B_h = B.^h;

F = zeros(N,N); 
% F is the matrix of the functional response for i eating j
for j = 1:N
    for i = 1:S % floral rewards don't eat
        if I(i,j) == 1
            F(i,j) = (omega(i,j) * B_h(j))/(((I(i,:) .* omega(i,:)) * B_h) + B0_h(j));
        end
    end
end

% plant competition
competition = calc_competition(competition_model,B,plants,K); 

consumption_rate = zeros(N,N);
% consumption_rate is the matrix of the rate at which i consumes j
for j = 1:N
    for i = 1:S % floral rewards don't eat
        if I(i,j) == 1
            consumption_rate(i,j) = x(i) * yij(i,j) * B(i) * F(i,j);
        end
    end
end

% functional response for accrual of benefits from pollination (assume
realized_benefit = calc_pollination_benefit(benefit_model,consumption_rate,pollinators,rewards,benefit_coefficient);

% pre-allocate flows for each species
production = zeros(N,1); % production is total growth minus metabolic losses
consumption_on = zeros(N,1); % how much loss due to consumption each guild is sustaining *not traditional
consumption_by = zeros(N,1); % how much consumption each guild is doing
metabolic_loss = zeros(N,1); % metabolic loss (0 for plants w/o pollinators) 

% ASSEMBLE THE FLOWS FOR EACH SPECIES
for i = 1:numel(rewards)
    production(rewards(i)) = beta * B(app(i)) - s * B(rewards(i));
    consumption_on(rewards(i)) = eij_eating_rewards * sum(consumption_rate(:,rewards(i))); % 
end

for i = 1:numel(app)
    production(app(i)) = competition * r_app * B(app(i)) * realized_benefit(i) - kappa * (beta * B(app(i)) - s * B(rewards(i)));
    consumption_on(app(i)) = eij_eating_vegetation *  sum(consumption_rate(:,app(i))); % 
    metabolic_loss(app(i)) = kappa * (beta * B(app(i)) - s * B(rewards(i)));
end

for i = 1:numel(wind)
    production(wind(i)) = competition * r_wind * B(wind(i));
    consumption_on(wind(i)) = eij_eating_vegetation * sum(consumption_rate(:,wind(i))); % 
end

for i = 1:numel(consumers)
    production(consumers(i)) = sum(consumption_rate(consumers(i),:)) - x(consumers(i)) * B(consumers(i));
    consumption_on(consumers(i)) = eij_eating_animals * sum(consumption_rate(:,consumers(i))); % 
    consumption_by(consumers(i)) = eij_eating_animals * sum(consumption_rate(consumers(i),consumers))...
        + eij_eating_vegetation * sum(consumption_rate(consumers(i),plants)) ...
        + eij_eating_rewards * sum(consumption_rate(consumers(i),rewards));
    metabolic_loss(consumers(i)) = x(consumers(i)) * B(consumers(i));
end

end