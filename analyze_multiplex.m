% This script analyzes outputs of Rewards Only Multiplex simulations.
% Specifically, it calculates DIVERSITY, PERSISTENCE, ABUNDANCE, and final
% FLOWS (including consumption, production, and metabolic loss) for the
% community and each guild of species at the end of simulations. 
% It generates a single file: Results_RO_Flows.txt with column headings 
% describing each entry and 24,276 rows of results for each network.
% Comments in the script below clarify what each entry is. 
% If an output file is not in the current directory, this script adds its 
% number to a list called "not_run" and removes those simulations from the 
% results table before saving. 
% Find & replace "RO" with "RP" for the Rewards Plus FW treatment.
%
% This code should take less than 30 mins to complete. 
% 
% CITE THIS CODE AS FOLLOWS:
% Hale, K.R.S. (2020). Pollinators in food webs: Mutualistic interactions 
%   increase diversity, stability, and function in multiplex networks

disp('Loading network structures and parameters...')

% Rewards Only (RO) Multiplex
load('networks_RO.mat')
load('metabolics_RO.mat')
load('simulation_parameters.mat','parameter_set')

networks = networks_RO;
metabolics = metabolics_RO;
parameter_set = parameter_set.multiplex;

num_networks = numel(networks);
num_samples = numel(parameter_set);
num_simulations = num_networks * num_samples;

% pre-allocate for results
simID = zeros(num_simulations,1);
S = zeros(num_simulations,1);
beta = zeros(num_simulations,1); 
diversity = zeros(num_simulations,1);
persistence = zeros(num_simulations,1);
community_biomass = zeros(num_simulations,1);

guild_results = zeros(num_simulations,70);

not_run = [];

for i = 1:num_simulations

    % load in outputs for analysis
    try
        filename = strcat('output',sprintf('%05d',i),'.mat');
        load(filename)
    
    disp(['Analyzing simulation: ', sprintf('%05d',i)])
    
    % match up network structure, metabolics, and parameters to outputs
    this_network = mod(i,num_networks);
    this_sample = floor(i/(num_networks)) + 1;
    
    if (this_network == 0)
        this_network = num_networks;
        this_sample = this_sample - 1;
        disp(i)
    end
    
    these_metabolics = metabolics(this_network);
    this_network = networks(this_network);
    this_sample = parameter_set(this_sample);
    
    % treatment
    simID(i,1) = i;
    if isempty(this_sample.beta)
        beta(i,1) = 0;
    else
        beta(i,1) = this_sample.beta;
    end
        
    S(i,1) = this_network.S; % number of species
    extinction_threshold = this_sample.extinction_threshold;

    % begin analysis
    diversity(i,1) = S(i,1) * output.persistence; 
    persistence(i,1) = output.persistence;
    
    final_biomasses = output.final_biomasses;
    community_biomass(i,1) = sum(final_biomasses,'omitnan'); 
    
    % calculate flows, evaluated at the final biomasses
    [total_production, total_consumption_on,...
        total_consumption_by, total_metabolic_loss, realized_benefit] = ...
        calc_multiplex_flows(this_network, these_metabolics, this_sample,...
        final_biomasses);
    
    % begin guild-specific analysis
    % declare guilds
    guilds(1).guild = this_network.wind; % plants w/o pollinators (e.g. wind-pollinated, selfing, etc.)
    guilds(2).guild = this_network.app; % plants w/ pollinators vegetation
    guilds(3).guild = this_network.rewards; % animal-pollinated rewards
    guilds(4).guild = this_network.plants; % all vegetation
    guilds(5).guild = setdiff(this_network.herbivores,51:S(i)); % strict herbivores
    guilds(6).guild = 51:S(i); % added-TL2 (pollinators in multiplex)
    guilds(7).guild = setdiff(union(this_network.omnivores,this_network.mammal_herbs),guilds(6).guild); % non-pollinator herbivorous omnivores
    guilds(8).guild = this_network.carnivores; % carnivores

    % omnivorous pollinators
    guilds(9).guild = intersect(union(this_network.omnivores,this_network.mammal_herbs),guilds(6).guild);
    % strictly herbivorous pollinators
    guilds(10).guild = setdiff(guilds(6).guild,guilds(9).guild);
    
    for j = 1:10
        
        % initial guild diversity
        i_guild_diversity = numel(guilds(j).guild);
        
        % guild diversity
        guild_diversity = sum(final_biomasses(guilds(j).guild) > extinction_threshold);
        
        % guild persistence
        guild_persistence = sum(final_biomasses(guilds(j).guild) > extinction_threshold)/numel(guilds(j).guild);
        
        % biomass distribution 
        guild_biomass = sum(final_biomasses(guilds(j).guild),'omitnan');
        
        % average guild flows
        guild_consumption_by = sum(total_consumption_by(guilds(j).guild),'omitnan');
        guild_metabolic_loss = sum(total_metabolic_loss(guilds(j).guild),'omitnan');
        guild_productivity = sum(total_production(guilds(j).guild),'omitnan');

        % log results
        guild_results(i,j) = i_guild_diversity;
        guild_results(i,j+10) = guild_diversity;
        guild_results(i,j+20) = guild_persistence;
        guild_results(i,j+30) = guild_biomass;
        guild_results(i,j+40) = guild_productivity;
        guild_results(i,j+50) = guild_consumption_by;
        guild_results(i,j+60) = guild_metabolic_loss;
        
    end
    
    % otherwise, add them to the not_run list
    catch
        not_run = [not_run; i];
        disp(i)
    end
    
end

% format results as a table
results = [simID, beta, S, diversity, persistence, community_biomass, guild_results];

% remove simulations that weren't run
results(not_run,:) = [];

% set up labels
VariableNames = {'simID','beta','S','final_diversity','persistence','community_biomass',...
    'i_div_wind','i_div_app_veg','i_div_app_rewards','i_div_all_vegetation','i_div_strict_herbs','i_div_addTL2','i_div_omni','i_div_carn','i_div_omni_poll','i_div_herb_poll',... % initial diversity
    'div_wind','div_app_veg','div_app_rewards','div_all_vegetation','div_strict_herbs','div_addTL2','div_omni','div_carn','div_omni_poll','div_herb_poll',... % final diversity
    'persist_wind','persist_app_veg','persist_app_rewards','persist_all_vegetation','persist_strict_herbs','persist_addTL2','persist_omni','persist_carn','persist_omni_poll','persist_herb_poll',... % persistence
    'bio_wind','bio_app_veg','bio_app_rewards','bio_all_vegetation','bio_strict_herbs','bio_addTL2','bio_omni','bio_carn','bio_omni_poll','bio_herb_poll',... % abundance (biomass)
    'prod_wind','prod_app_veg','prod_app_rewards','prod_all_vegetation','prod_strict_herbs','prod_addTL2','prod_omni','prod_carn','prod_omni_poll','prod_herb_poll',... % productivity 
    'cons_by_wind','cons_by_app_veg','cons_by_app_rewards','cons_by_all_vegetation','cons_by_strict_herbs','cons_by_addTL2','cons_by_omni','cons_by_carn','cons_by_omni_poll','cons_by_herb_poll',... % consumption rate
    'loss_wind','loss_app_veg','loss_app_rewards','loss_all_vegetation','loss_strict_herbs','loss_addTL2','loss_omni','loss_carn','loss_omni_poll','loss_herb_poll'}; % metabolic loss
    
results_table = array2table(results,'VariableNames',VariableNames);

% save results
writetable(results_table,'Results_RO_Flows.txt')