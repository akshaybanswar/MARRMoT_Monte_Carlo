function [ a,b,c,d,e ] = MonteCarloSamplingResults( numSample,model,input_climate,input_solver )
%MonteCarloSamplingResults runs a Monte Carlo Analysis of a model and generates the
%sampling plots.
%
% Input:
% numSample         - the requested number of samples
% model             - model name
% input_climate     - time series of climatic fluxes in simulation period
% input_solver      - solver settings
%
% Output:
% Returning the Output
% a = input_theta;
% b = output_ex;
% c = output_in;
% d = output_ss;
% e = output_waterbalance;

% Display progress
% -------------------------------------------------------------------------
disp('Monte Carlo Analysis Started')
disp(' ')

% Prepare the model
% -------------------------------------------------------------------------

% Extract the parameter range for this model
% NOTE: files with parameter ranges are found in the folder
% './MARRMoT/Models/Parameter range files/'
model_range = feval([model,'_parameter_ranges']);                          % Call the function that stores parameter ranges

% Find number of stores and parameters
numPar    = size(model_range,1);
numStore    = str2double(model(end-1));

% Set the inital storages
% NOTE: see the model function for the order in which stores are given. For
% HyMOD, this is on lines 86-91.
input_s0    = zeros(numStore,1);

% GENERATE SAMPLE
input_theta = (model_range(:,1)+lhsdesign(numPar,numSample).*(model_range(:,2)-model_range(:,1)))';

% Pre-allocating NaN
for j = 1:numSample
    output_ex{1,1}{1,j}.Q(1:length(input_climate.date)) = NaN;
    output_ex{1,1}{1,j}.Ea(1:length(input_climate.date)) = NaN;
end

% Monte Carlo Simulation
parfor i = 1:numSample

    % Run the model
    [output_ex{i},...                                                      % Fluxes leaving the model: simulated flow (Q) and evaporation (Ea)
        output_in{i},...                                                   % Internal model fluxes
        output_ss{i},....                                                  % Internal storages
        output_waterbalance{i}] = ...                                      % Water balance check
        feval(model,...                                                    % Model function name
        input_climate,...                                                  % Time series of climatic fluxes in simulation period
        input_s0,...                                                       % Initial storages
        input_theta(i,:),...                                               % input_theta{i},...   % Parameter values
        input_solver);                                                     % Details of numerical time-stepping scheme

    disp(sprintf('%s iteration number : %d',model,i));

end

% Returning the Output
a = input_theta;
b = output_ex;
c = output_in;
d = output_ss;
e = output_waterbalance;

% Display progress
% -------------------------------------------------------------------------
disp(' ')
disp('Monte Carlo Analysis Finished')
end

