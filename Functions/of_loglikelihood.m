function [LL] = of_loglikelihood(obs,sim,warmup,varargin)
%% check inputs and set defaults
if nargin < 3
    error('Not enugh input arguments')
elseif nargin > 4
    error('Too many inputs.')

elseif nargin == 3
    w = [1,1,1];                                                           % no weights specified, use defaults

elseif nargin == 4
    if size(varargin{1}) == [1,3] | size(varargin{1}) == [3,1]             % check weights variable for size
        w = varargin{1};                                                   % apply weights if size = correct
    else
        error('Weights should be a 3x1 or 1x3 vector.')                    % or throw error
    end
end

% check time series size and rotate one if needed
if checkTimeseriesSize(obs,sim) == 0
    error('Time series not of equal size.')
elseif checkTimeseriesSize(obs,sim) == 2
    sim = sim';                                                             % 2 indicates that obs and sim are the same size but have different orientations
end

% ##################
obs = obs(warmup:end);
sim = sim(warmup:end);
% ##################

%% check for missing values
% -999 is used to denote missing values in observed data, but this is later
% scaled by area. Therefore we check for all negative values, and ignore those.
idx = find(obs >= 0);

%% Calculation of AIC & BIC by Residual Sum of Squared Error
diff = obs(idx)-sim(idx);
diff_square = diff.^2;
RSS = sum(diff_square);
numObs = length(idx);
meanobs = mean(obs(idx));
DiffFromMean = obs(idx)-meanobs;
DiffFromMean_square = DiffFromMean.^2;
ErrorVariance = (sum(DiffFromMean_square))/numObs;
LL = log(RSS/ErrorVariance);
%LL = -((numObs/2)*log(2*pi))-((numObs/2)*log(var(obs(idx))))-(((var(obs(idx)))/2))*(diff_square);
end