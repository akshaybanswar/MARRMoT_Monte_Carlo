function [bic] = of_bic(obs,sim,warmup,numPar,varargin)
%% check inputs and set defaults
if nargin < 4
    error('Not enugh input arguments')
elseif nargin > 5
    error('Too many inputs.')

elseif nargin == 4
    w = [1,1,1];                                                           % no weights specified, use defaults

elseif nargin == 5
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

bic = numObs * log(RSS/numObs) + numPar * log(numObs);

end