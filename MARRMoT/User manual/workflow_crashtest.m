% This file is part of the Modular Assessment of Rainfall-Runoff Models 
% Toolbox (MARRMoT) � User manual. It contains an example of crash tests 
% that can be used to quality control a newly create model function. See 
% section 4 in the manual for details.
% 
% Author:   Wouter J.M. Knoben
% Date:     29-09-2018
% Contact:  w.j.m.knoben@bristol.ac.uk
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.

% This example workflow includes 4 steps:
%
% 1. Specify the model and solver
% 2. Load test data
% 3. Parameter extremes crasht test
% 4. Random parameter value water balance check

%% 1. Specify model and prepare some auxiliary variables
% Select model
model  = 'm_nn_example_7p_3s';                     % Model function

% Extract parameters from model parameter function
theta  = feval([model,'_parameter_ranges']);   

% Find number of stores, parameters and define initial storages
numPar = size(theta,1);
numSto = str2double(model(end-1));
sIni   = zeros(numSto,1);

% Specify solver settings
solver.name              = 'createOdeApprox_IE';   % Use Implicit Euler to approximate ODE's
solver.resnorm_tolerance = 0.1;                    % Convergence tolerance
solver.resnorm_maxiter   = 5;                      % Maximum number of re-runs

%% 2. Prepare test data
% Load data
load MARRMoT_example_data.mat

% Create a climatology data input structure. 
% Using 1 year speeds up the process 
input_climatology.precip  = data_MARRMoT_examples.precipitation(1:365); 
input_climatology.temp    = data_MARRMoT_examples.temperature(1:365);
input_climatology.pet     = data_MARRMoT_examples.potential_evapotranspiration(1:365);
input_climatology.delta_t = 1;                                              % [d]

%% 3. Run a parameter extreme crash test
% Find all possible permutations of min/max parameter values
combos = dec2bin(0:(2^numPar-1))-'0'+1;

% Create an empty vector to store water balance checks
check_wb = NaN.*zeros(length(combos),1);

% Crash test the model
for i = 1:length(combos)
    
    % Display a counter
    disp(['Currently on ',num2str(i),'/',num2str(length(combos))]);
    
    % Get parameters
    theta_select = theta(sub2ind(size(theta),[1:numPar],combos(i,:)));
    
    % Run model
    % The first three outputs are structures with flux and store time
    % series. These are not required now. The fourth output is the water
    % balance variable. Saving and plotting this helps with models that
    % have many parameters, because it is infeasible to follow the output
    % to the screen for larger models.
    [~,~,~,check_wb(i)] = feval(model,...               % Model name
                                input_climatology,...   % Time series of climatic fluxes in simulation period
                                sIni,...                % Initial guess of store values
                                theta_select,...        % Parameter values
                                solver);                % Numerical scheme to be used
                   
end

% Check the water balance errors
figure('color','w')
    plot(check_wb)
    xlabel('Parameter set')
    ylabel('Water balance error [mm]')
   
%% 4. Random parameter value water balance check
% Specify the number of samples
n = 50;

% Create an empty vector to store water balance checks 
check_wb = NaN.*zeros(n,1);

% Create an empty matrix to store randomized parameter sets, in case
% further investigation is needed
theta_select = NaN.*zeros(n,numPar);

for i = 1:n

    % Display a counter
    disp(['Currently on ',num2str(i),'/',num2str(n)]);
    
    % Get parameters
    theta_select(i,:) = theta(:,1)+rand(numPar,1).*(theta(:,2)-theta(:,1));
    
    % Run model
    [~,~,~,check_wb(i)] = feval(model,...               % Model name
                                input_climatology,...   % Time series of climatic fluxes in simulation period
                                sIni,...                % Initial guess of store values
                                theta_select(i,:),...   % Parameter values
                                solver);                % Numerical scheme to be used
end

% Check the water balance errors
figure('color','w')
    plot(check_wb)
    xlabel('Parameter set')
    ylabel('Water balance error [mm]')
