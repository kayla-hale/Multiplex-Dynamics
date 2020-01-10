function [omega] = calc_preference(preference_model,I)
%CALC_PREFERENCE Calculate the preference each consumer has for each of its
%resources. Two models are implemented, corresponding to the generalist 
%models described in Williams 2008: 'weak', in which each consumer splits 
%its foraging time equally among all its resources even if they go 
%extinct, and 'strong' where each consumer has equal preference for all 
%its resources so that it forages proportionally to their abundance. Let N
%be the length of IR. 
%   omega is an NxN matrix of doubles corresponding to row i's preference
%       for the resource in column j
%   preference_model is a string specifying how preference is modeled:
%       'weak' or 'strong'
%   I is an NxN logical matrix where 1 indicates that row i interacts with
%       column j and 0 indicates no interaction
%
% CITE THIS CODE AS FOLLOWS:
% Hale, K.R.S. (2020). Pollinators in food webs?Mutualistic interactions 
%   increase diversity, stability, and function in multiplex networks

    switch(preference_model)
        case 'weak'
            % total number of interaction partners for species in row i
            diet_sizes = sum(I,2);
            omega = 1./diet_sizes;
            omega(isinf(omega)) = 0;
            omega = I .* omega; 
        
        case 'strong'
            omega = I;
            
        otherwise
            error('Requested preference model not implemented!')
    end
        
end

