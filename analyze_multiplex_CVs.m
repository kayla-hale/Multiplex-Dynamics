% This script analyzes solutions of Rewards Only Multiplex simulations.
% Specifically, it calculates BIOMASS VARIABILITY for the final community,
% each guild, and species at the end of simulations. Vegetation and Rewards
% of plants w/ pollinators are calculated both as separate guild and are
% summed into one guild. 
% It generates a single file: Results_RO_CVs.txt with column headings 
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

% Rewards Only (RO) Multiplex Networks
load('networks_RO.mat')
load('metabolics_RO.mat')
load('simulation_parameters.mat')

networks = networks_RO;
metabolics = metabolics_RO;
parameter_set = parameter_set.multiplex;

num_networks = numel(networks);
num_samples = numel(parameter_set);
num_simulations = num_networks * num_samples;
num_guilds = 10;

% pre-allocate for results
simID = zeros(num_simulations,1);
S = zeros(num_simulations,1);
beta = zeros(num_simulations,1);
persistence = zeros(num_simulations,1);

avg_species_CV = zeros(num_simulations,1);
avg_guild_CV_stzd = zeros(num_simulations,1);
community_CV = zeros(num_simulations,1);

guild_results = zeros(num_simulations,(num_guilds+1)*2);

not_run = [];

for i = 400:num_simulations
    
    % load in outputs for analysis
    try
        filename = strcat('solution',sprintf('%05d',i),'.mat');
        load(filename)
        
    disp(['Analyzing simulation: ', sprintf('%05d',i)])
    
    % match up network structure, metabolics, and parameters to outputs
    this_network = mod(i,num_networks);
    this_sample = floor(i/(num_networks)) + 1;
    
    if (this_network == 0)
        this_network = num_networks;
        this_sample = this_sample - 1;
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
    this_S = this_network.S;
    extinction_threshold = this_sample.extinction_threshold;

    % begin analysis
    % species persistence
    surviving_species = solution(1:this_S,end) > extinction_threshold;
    persistence(i,1) = sum(surviving_species)/this_S;

    solution = solution(:,end-1000:end);

    solution(~surviving_species,:) = NaN; 

    % coefficients of variation (of species that aren't extinct throughout the interval) 
    % consider animal-pollinated plants' vegetation and rewards together
    % (to standardize division by surviving S between multiplex & FWs)
    solution_S = solution(1:this_S,:);
    solution_S(this_network.app,:) = solution(this_network.rewards,:) + solution(this_network.app,:);
    % each species' CV
    each_species = std(solution_S,[],2)./mean(solution_S,2); % extinct species are NaN
    % average (surviving) species' CV
    avg_species_CV(i,1) = mean(each_species,'omitnan'); % omit extinct species
    % community CV (total summed biomass)
    community_time_series = sum(solution_S,1,'omitnan'); % omit extinct species
    community_CV(i,1) = std(community_time_series)/mean(community_time_series); % extinct already omitted

    % begin guild-specific analysis
    % declare guilds
    guilds(1).guild = this_network.wind; % wind-pollinated
    guilds(2).guild = this_network.app; % animal-pollinated vegetation
    guilds(3).guild = this_network.rewards; % animal-pollinated rewards
    guilds(4).guild = this_network.plants; % all vegetation
    guilds(5).guild = setdiff(this_network.herbivores,51:this_S); % strict herbivores
    guilds(6).guild = 51:this_S; % added-TL2 (pollinators in multiplex)
    guilds(7).guild = setdiff(union(this_network.omnivores,this_network.mammal_herbs),guilds(6).guild); % non-pollinator herbivorous omnivores
    guilds(8).guild = this_network.carnivores; % carnivores

    % omnivorous pollinators
    guilds(9).guild = intersect(union(this_network.omnivores,this_network.mammal_herbs),guilds(6).guild);
    % strictly herbivorous pollinators
    guilds(10).guild = setdiff(guilds(6).guild,guilds(9).guild);

    for k = 1:num_guilds

        each_guild_species_time_series = solution(guilds(k).guild,:);

        % average within guild (surviving) species' CV
        each_guild_species_CV = std(each_guild_species_time_series,[],2)./mean(each_guild_species_time_series,2); % NaN's are extinct species
        avg_guild_species_CV = mean(each_guild_species_CV,'omitnan'); % omit extinct species

        % guild CV (total summed biomass)
        guild_time_series = sum(each_guild_species_time_series,1,'omitnan'); % omit extinct species
        guild_CV = std(guild_time_series)/mean(guild_time_series); % NaNs are already omitted

        % log results
        guild_results(i,k) = avg_guild_species_CV; 
        guild_results(i,k+num_guilds+1) = guild_CV;

    end

    % animal-pollinated vegetation & rewards together (stored in
    % solution_S(this_network.app,:))

    % average within guild (surviving) species' CV
    each_guild_species_time_series = solution_S(this_network.app,:);
    each_guild_species_CV = std(each_guild_species_time_series,[],2)./mean(each_guild_species_time_series,2); % NaN's are extinct species
    avg_guild_species_CV = mean(each_guild_species_CV,'omitnan'); % omit extinct species

    % guild CV (total summed biomass)
    guild_time_series = sum(each_guild_species_time_series,1,'omitnan'); % omit extinct species
    guild_CV = std(guild_time_series)/mean(guild_time_series); % NaNs are already omitted

    guild_results(i,(num_guilds+1)) = avg_guild_species_CV; 
    guild_results(i,2*(num_guilds+1)) = guild_CV;


    % average guild CV, standardized by max. number of species in each 
    % guild for comparison between multiplex and FW simulations
    avg_guild_CV_stzd(i,1) = sum(guild_results(i,4+num_guilds+1:8+num_guilds+1),'omitnan')/5;

    % otherwise, add them to the not_run list
    catch
        disp(i)
        not_run = [not_run; i];
    end

end

% format results as a table
results = [simID, beta, S, persistence, avg_species_CV, avg_guild_CV_stzd, community_CV, guild_results];

% remove simulations that weren't run
results(not_run,:) = [];

% set up labels
VariableNames = {'simID','beta','S','persistence','avg_species_CV','avg_guild_CV_stzd','community_CV',...
    'wind_species_CV','app_veg_species_CV','app_rewards_species_CV','all_vegetation_species_CV','strict_herbs_species_CV','addTL2_species_CV','omni_species_CV','carn_species_CV','omni_poll_species_CV','herb_poll_species_CV','app_all_species_CV',...
    'wind_guild_CV','app_veg_guild_CV','app_rewards_guild_CV','all_vegetation_guild_CV','strict_herbs_guild_CV','addTL2_guild_CV','omni_guild_CV','carn_guild_CV','omni_poll_guild_CV','herb_poll_guild_CV','app_all_guild_CV'};

results_table = array2table(results,'VariableNames',VariableNames);

% save simulation summary table
writetable(results_table,'Summary_RO_CVs.txt')
