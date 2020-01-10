# Multiplex-Dynamics
Companion code for Hale et al. 2020: Pollinators in food webs: Mutualistic interactions increase diversity, stability, and function in multiplex networks

% All code in this folder is companion to the manuscript: 
% Hale et al. 2020 Pollinators in food webs: Mutualistic interactions 
%   increase diversity, stability, and function in multiplex networks
% See the main-text Methods for a detailed description of the dynamic 
%   models and analysis implemented in this code.
%
% SYSTEM REQUIREMENTS
% Hardware Requirements
% 	A standard laptop computer should suffice to run this code. It has
%	been tested with a MacBook with 2.2 GHz Intel i7 Processor and 16 
%	GB Memory.
%
% Software Requirements
%	This code is supported for MATLAB R2018b or later. It has been 
%	tested on macOS Mojave 10.14.6. Some functions use the Parallel
%   	Computing Toolbox, but that is optional, as described by comments 
%   	in the scripts.
%
%
% INSTALLATION
% 	To use this code, download the whole folder and use it as a
% 	working directory for MATLAB R2018b or later. Execute
% 	RUN_DEMO_SIMULATIONS.m. The code is working correctly if it produces
% 	six figures that match the panels in Fig. S1 of the manuscript 
%	(expect approximately 2 minutes for run time). 
%
% REPRODUCTION
%	To fully reproduce the main-text results, use the
%	RUN_MULTIPLEX_SIMULATIONS.m & RUN_FOOD_WEB_SIMULATIONS.m scripts
% 	as specified in their description (i.e. for both "RO" and "RP"
%	treatments). Then, use ANALYZE_MULTIPLEX.m & ANALYZE_MULTIPLEX_CVS.m
%	scripts for the multiplex simulations or ANALYZE_FOOD_WEB.m & 
%	ANALYZE_FOOD_WEB_CVS.m to generate space-separated .txt files of all
% 	outputs used in the main-text figures. These functions run & analyze
%	dynamic simulations given network structures and parameter sets as 
%	inputs, which are included in this folder. See below for description.
%	The Methods section and Figure legends describe how summary statistics
%	including means and standard deviations across treatment and for each
%	level of initial diversity within treatments were calculated. These
%	calculations are easily performed in JMP (version Pro 14 was used to
%	create main-text Figs. 4 & 5).
%
% LICENSE
% 	The MIT License (MIT)
% 	Copyright (c) 2020 Kayla R. S. Hale
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
% copies of the Software, and to permit persons to whom the Software is 
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
% SOFTWARE.
%
% CITE THIS CODE AS FOLLOWS:
% Hale, K.R.S. (2020). Pollinators in food webs: Mutualistic interactions 
%   increase diversity, stability, and function in multiplex networks
%
% Other Notes
% Individual functions are quite thoroughly commented. The remainder of 
% 	this README describes the inter-dependencies between the functions and 
%	the nature of the data/input files. 
%

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% RUN_DEMO_SIMULATIONS.m takes a single network subjected to all six 
% main-text treatments (RO FW, RO Low, RO High, RP FW, RP Low, RP High), 
% simulates dynamics, and displays the associated timeseries and persistence.
% These simulations correspond to those shown in Fig. S1. I
% recommend exploring these timeseries and the code that generated them to
% familiarize yourself with the model used in this paper. A description of
% the functions and other inputs RUN_DEMO_SIMULATIONS needs is given below.
% This demo script and the other RUN scripts load in the needed data and 
% call the needed functions automatically when the files included in this 
% folder are in your working directory. 

% An initial network structure and a set of parameters is needed to execute 
% the dynamic models. The 24,276 network structures for each treatment
% are stored in the NETWORKS .mat files. The files are a 24276x1 array of 
% structs. Each struct has 27 fields that describe the structural properties
% of the network (not all of them necessary for this manuscript).
%   fileID: a string describing the pollination sub-network (p) and niche  
%       model (f) used to construct the whole network (see Main-Text Fig. 2 
%       and Methods S1 for details of network construction).
%   S: number of species (initial diversity)
%   N: number of state variables including rewards
%   feeding_C: connectance of the niche model food web
%   pollination_C: connectance of the pollination sub-network
%   num_plants: number of plants
%   fraction_app: fraction of the plants that have pollinators (are 
%       animal-pollinated plants, hereafter referred to as "app")
%   pnest: value used in ThÈbault & Fontaine (2010)'s algorithm to generate
%       pollination networks (see Methods S2)
%   NODF: nestedness of pollination network
%   degree_heterogeneity: degree heterogeneity of pollination network
%   fraction_omnivores: fraction of consumer species that feed on both 
%       plants and animals in the initial niche model
%   n: niche values from the niche model
%   I: NxN binary interaction network; species are in entries 1:S, rewards
%       nodes are in S+1:N. Food web networks are SxS. 
%   plants, app, rewards, etc.: vectors of indices identifying species in 
%       that guild. If a guild is empty, the vector is empty also (== []).
%	In particular, TL2 is the union of herbivores and added-animals.
%   pollinator_list: an array of structs, with each struct a list of
%       pollinators that link to that app's rewards node
%   C: connectance of the entire network (see Fig. S4)
%
% The METABOLICS .mat files provide metabolic rate, which is used to 
% calculate allometrically-scaled parameters in the DYNAMICS functions. The 
% files are a 24276x1 array of structs, with each struct corresponding to 
% the network with the same fileID in the corresponding NETWORKS file. Each 
% struct has 6 fields.
%   fileID: same as in the NETWORKS files
%   S: number of state variables including rewards (i.e. corresponds to N
%       in the NETWORKS files)
%   mets_algorithm: string used to cite the algorithm used to construct 
%       metabolic rates: 'AAAI2012' refers to Martinez et al. 2012, published
%       in the journal AAAI; also see Methods & Table S4 for details
%   swTL: Sx1 double vector of short-weighted trophic levels for each node;
%       plants & rewards are defined to have swTL == 1
%   x: metabolic rate for each node; plants & rewards are defined to have x
%       == 1
%   web_strategy: string used to identify consumers' strategy for calculat-
%       ing allometric rates: only 'invertebrate' was used in this study
%
% Finally, the SIMULATION_PARAMETERS.mat struct specifies all other (non-
% allometric) parameters. See comments in DYNAMICS functions or Table S4
% for details. SIMULATION_PARAMETERS loads in a struct with field 
% .multiplex with the two (Low & High rewards productivity) multiplex 
% parameterizations and .food_web with the single FW parameterization. 

% To simulate the dynamics of Food Webs (FW), the function
% FOOD_WEB_DYNAMICS.m takes an initial network structure and a set of
% parameters including those determined by allometry. It then numerically
% integrates a set of differential equations (see Methods) over 5000
% timesteps using ode45. It uses the simple functions calc_competition.m 
% and calc_preference.m during execution. It saves two files: 
% OUTPUT files are  .mat files containing persistence, details of the
% treatment and network, biomass distribution at the end of the simulation,
% and when species went extinct if they did. 
% SOLUTION files are .mat files containing biomass at each timestep for 
% each species.

% To simulate the dynamics of Multiplex networks, the function
% MULTIPLEX_DYNAMICS.m takes an initial network structure and a set of
% parameters including those determined by allometry. It then numerically
% integrates a set of differential equations (see Methods) over 5000
% timesteps using ode15s. In addition to calc_competition.m and 
% calc_preference.m, it also uses calc_pollination_benefit during execution.
% It also saves OUTPUT and SOLUTION files. 

% The other RUN_ scripts (RUN_FOOD_WEB_SIMULATIONS.m and 
% RUN_MULTIPLEX_SIMULATIONS.m) are the scripts used to generate all results 
% in the manuscript. They are written so that multiple simulations can be 
% run in parallel which requires the Parallel Computing Toolbox. This is 
% optional, but speeds up execution time substantially. All RUN_ scripts
% call the same functions to simulate network dynamics. Simulation analyses 
% are performed by the ANALYZE scripts. All analyses except for species and 
% guild variability are performed in the ANALYZE_FOOD_WEB.m or
% ANALYZE_MULTIPLEX.m scripts by reading in all OUTPUT .mat files generated
% from the RUN_ scripts. Variability is calculated with the ANALYZE_ CVs.m
% scripts by reading in the SOLUTION .mat files. 
