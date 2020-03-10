function [competition] = calc_competition(competition_model,B,plants,K,sum_R)
%CALC_COMPETITION Calculate the effect of competition on plant per-capita
%growth rate. Two models are implemented: 'community-wide', where all
%plant species' vegetative biomass (B) contribute to  how close the plant
%community is to carrying capacity, and 'community-wide_control', which is
%used for rewards partitioning control simulations where plants' biomass 
%(B) includes rewards, so, to match the original simulations, the sum of 
%all species' rewards biomass from the original simulations is subtracted 
%off the competition calculation. Using sum_R = 0 or sum_R = [] makes  
%'community-wide_control' identical to the 'community-wide' calculation. 
%
%   competition is a double between 0 and 1 that is equal to the realized
%       effect of competition all plants in the community exert on each
%       species
%   competition_model is a string specifying how plant competition is
%       modeled: 'commmunity-wide' or 'community-wide_control'
%   B is a vector of doubles equal to biomasses of each species
%   plants is an integer vector of indices (corresponding to which entries 
%       are plants in B)
%   K is a double equal to the community-wide carrying capacity of plants
%   sum_R is a double equal to summed rewards over all plants (averaged
%       over time from the original simulations)
%
% CITE THIS CODE AS FOLLOWS:
% Hale, K.R.S. (2020). Mutualistic interactions increase diversity, 
%   stability, and function in multiplex networks of pollinators in food webs

    switch(competition_model)
        case 'community-wide'
            competition = 1 - sum(B(plants))/K;
            
        case 'community-wide_control'
            competition = 1 - sum(B(plants))/K + sum_R/K;
            
        otherwise
            error('Requested competition model not implemented!')
    end

end