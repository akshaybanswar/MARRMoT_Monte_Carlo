function [ out,UH ] = uh_3_half( in, d_base, delta_t )
%uh_3_half Unit Hydrograph [days] with half a triangle (linear)
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
%   Inputs
%   in      - volume to be routed
%   d_base  - time base of routing delay [d]
%   delta_t - time step size [d]    
%
%   Unit hydrograph spreads the input volume over a time period delay.
%   Percentage of input returned only increases. 
%   I.e. d_base = 3.8 [days], delta_t = 1:
%   UH(1) = 0.04  [% of inflow]
%   UH(2) = 0.17
%   UH(3) = 0.35
%   UH(4) = 0.45

%%INPUTS
if any(size(in)) > 1; error('UH input should be a single value.'); end

%%TIME STEP SIZE
delay = d_base/delta_t;
if delay == 0; delay = 1; end   % any value below t = 1 means no delay, 
                                % but zero leads to problems
tt = 1:ceil(delay);             % time series for which we need UH 
                                % ordinates [days]

%%UNIT HYDROGRAPH
% The area under the unit hydrograph by definition sums to 1. Thus the area
% is S(t=0 to t = delay) t*[ff: fraction of flow to move per time step] dt
% Analytical solution is [1/2 * t^2 + c]*ff, with c = 0. Thus the fraction
% of flow step size is:
ff = 1/(0.5*delay^2);           

%%EMPTIES
UH = zeros(1,length(tt));

%%UNIT HYDROGRAPH
for t = 1:length(tt)
    if t <= delay
        UH(t) = ff.*(0.5*t^2 - 0.5*(t-1)^2);
    else
        UH(t) = ff.*(0.5*delay^2 - 0.5*(t-1)^2);
    end
end

%%DISPERSE VOLUME
out = in.*UH;

end

