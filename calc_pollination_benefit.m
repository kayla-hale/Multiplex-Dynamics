function [realized_pollination_benefit] = calc_pollination_benefit(benefit_model,consumption_rate,pollinators,rewards,benefit_coefficient)
%CALC_POLLINATION_BENEFIT Calculate the benefit due to pollination to the
%per-capita growth rate of animal-pollinated plants. Only one model is
%implemented so far: 'saturating', where the effect of reproductive
%services (accounting for the quality & quantity of pollination from each
%pollinator) saturate. 
%   realized_pollination_benefit is a vector of doubles bounded between 0 
%       and 1 equal to the amount of realized benefit to per-capita growth 
%       rate to each plant provided by all of its pollinators
%   competition_model is a string specifying how plant competition is
%       modeled: 'saturating'
%   consumption_rate is a double matrix where each entry corresponds to 
%       consumer in row i's consumption rate on the resource in column j
%   pollinators is an integer vector of indices (corresponding to which 
%       rows & columns are pollinators in consumption_rate)
%   rewards is an integer vector of indices (corresponding to which 
%       rows & columns are rewards in consumption_rate); rewards are
%       ordered matching the order of animal-pollinated plants' (hereafter
%       referred to as "app") indices
%   benefit_coefficient is a double that parameterizes the reproductive
%   	services functional response
%
% CITE THIS CODE AS FOLLOWS:
% Hale, K.R.S. (2020). Mutualistic interactions increase diversity, 
%   stability, and function in multiplex networks of pollinators in food webs
    
    switch(benefit_model)
        case 'saturating'
            reproductive_services_matrix = zeros(numel(rewards),numel(pollinators));

            for j = 1:numel(pollinators)
                for i = 1:numel(rewards)
                    % reprodutive services provided by pollinators(j) to
                    % app(i) is accounts for the quality and quantity of
                    % consumption by pollinators(j) on rewards(i) 
                    % i.e. reproductive services are a function of
                    % consumption rate by the pollinator on the plants'
                    % rewards
                    reproductive_services_matrix(i,j) = (consumption_rate(pollinators(j),rewards(i))^2 / sum(consumption_rate(pollinators(j),:)));
                end
            end
            
            reproductive_services_matrix(isnan(reproductive_services_matrix)) = 0;
            
            % total reproductive services provided to each plant
            reproductive_services = sum(reproductive_services_matrix,2);
            
            % realized_pollination_benefit: functional response for benefit 
            % accrual to app from pollinators (executed in dynamics)
            realized_pollination_benefit = reproductive_services./(benefit_coefficient + reproductive_services);
            
            if any(isnan(realized_pollination_benefit))
                error('NaNs in benefit function!/n')
            end
            
        otherwise
            error('Requested benefit model not implemented!/n')
    end
end


