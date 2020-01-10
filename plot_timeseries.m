function plot_timeseries(this_network,solution)
%PLOT_TIMESERIES Plot a log-log dynamical timeseries, with 
%species' abundances over time colored according to their guilds. The
%figure this function generates corresponds to the panels in Fig. 
%S1. 
%This function works for both multiplex
%networks (with rewards) and FWs (no rewards). It also prints the species
%persistence in the title of the figure, calculated assuming the extinction
%threshold is 10^(-6). The y-axis is set to a minimum of 10^(-4) & maximum
%of 500.
%   this_network is a structure including initial diversity (.S) & integer
%       indices of all guilds (e.g. .carnivores, .omnivores, etc.)
%   solution is an array where each column corresponds to a timestep
%       (starting at 0) and each row corresponds to a state variable (i.e.
%       species or rewards) whose biomass is changing over time
%
% CITE THIS CODE AS FOLLOWS:
% Hale, K.R.S. (2020). Pollinators in food webs?Mutualistic interactions 
%   increase diversity, stability, and function in multiplex networks

% initial diversity
S = this_network.S;
% number of state variables in the solution, = S + number of rewards nodes
N = size(solution,1);
% simulation length (number of timesteps)
x_max = size(solution,2);
x = (0:x_max-1)';
% species persistence
S_persist = sum(solution(1:S,end) > 10^(-6))/S;

% declare guilds
added_species = 51:S; % added-TL2 (pollinators in multiplex)
   
carnivores = this_network.carnivores; % carnivores
omnivores = setdiff(union(this_network.omnivores,this_network.mammal_herbs),added_species); % non-pollinator herbivorous omnivores
omni_polls = intersect(union(this_network.omnivores,this_network.mammal_herbs),added_species); % omnivorous pollinators
herb_polls = setdiff(added_species,omni_polls); % strictly herbivorous pollinators
herbivores = setdiff(this_network.herbivores,51:S); % strict herbivores
rewards = this_network.rewards; % animal-pollinated plants' rewards
app = this_network.app; % animal-pollinated plants' vegetation
wind = this_network.wind; % plants w/o pollinators (e.g. wind-pollinated, selfing, etc.)

% declare guilds' colors
colors = ...
    [102 0 0; ... % 1: carnivores
    255 0 0; ... % 2: omnivores
    255 153 0; ... % 3: omni polls
    204 255 51; ... % 4: herb polls
    51 255 204; ... % 5: herbivores
    153 51 204; ... % 6: rewards
    51 153 255; ... % 7: app
    0 0 255]./256; % 8: wind
    
figure

% to make the legend correctly
loglog(x(end), solution(wind(1),end),'LineWidth',2,'Color',colors(1,:))
hold on
loglog(x(end), solution(wind(1),end),'LineWidth',2,'Color',colors(2,:))
loglog(x(end), solution(wind(1),end),'LineWidth',2,'Color',colors(3,:))
loglog(x(end), solution(wind(1),end),'LineWidth',2,'Color',colors(4,:))
loglog(x(end), solution(wind(1),end),'LineWidth',2,'Color',colors(5,:))
loglog(x(end), solution(wind(1),end),'LineWidth',2,'Color',colors(6,:))
loglog(x(end), solution(wind(1),end),'LineWidth',2,'Color',colors(7,:))
loglog(x(end), solution(wind(1),end),'LineWidth',2,'Color',colors(8,:))

% plot the timeseries
if ~isempty(wind)
    loglog(x,solution(wind,:),'LineWidth',2,'Color',colors(8,:))
    hold on
end
if ~isempty(omnivores)
    loglog(x,solution(omnivores,:),'LineWidth',2,'Color',colors(2,:))
    hold on
end
if ~isempty(rewards) && N > S
    loglog(x,solution(rewards,:),'LineWidth',2,'Color',colors(6,:))
    hold on
end
if ~isempty(app)
    loglog(x,solution(app,:),'LineWidth',2,'Color',colors(7,:))
    hold on
end
if ~isempty(omni_polls)
    loglog(x,solution(omni_polls,:),'LineWidth',2,'Color',colors(3,:))
    hold on
end
if ~isempty(herb_polls)
    loglog(x,solution(herb_polls,:),'LineWidth',2,'Color',colors(4,:))
    hold on
end
if ~isempty(herbivores)
    loglog(x,solution(herbivores,:),'LineWidth',2,'Color',colors(5,:))
    hold on
end
if ~isempty(carnivores)
    loglog(x,solution(carnivores,:),'LineWidth',2,'Color',colors(1,:))
    hold on
end

xlabel('log10 Time Steps')
ylabel('log10 Biomass')
legend('Carnivores','Omnivores','+Omnivores/Pollinators','+Herbivores/Pollinators',...
    'Herbivores','Rewards','Plants w/ Pollinators','Plants w/o Pollinators'); 
grid on
% grid minor
title(['Persistence = ' num2str(S_persist*100,3) '%'])
xlim([1 x_max])
set(gcf,'color','w');
set(gca,'FontSize',16);
ylim([10^(-4) 500])
end
