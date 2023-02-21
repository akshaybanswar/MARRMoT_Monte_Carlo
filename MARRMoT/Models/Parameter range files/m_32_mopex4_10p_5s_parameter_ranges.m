function [ theta ] = m_32_mopex4_10p_5s_parameter_ranges( )
%m_32_mopex4_10p_5s_parameter_ranges Provides parameter ranges for calibration
%   of the 5-store MOPEX-4 model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% Ye, S., Yaeger, M., Coopersmith, E., Cheng, L., & Sivapalan, M. (2012). 
% Exploring the physical controls of regional patterns of flow duration 
% curves - Part 2: Role of seasonality, the regime curve, and associated 
% process controls. Hydrology and Earth System Sciences, 16(11), 4447�4465.
% http://doi.org/10.5194/hess-16-4447-2012

theta = [-3  , 3;       % tcrit, Snowfall & snowmelt temperature [oC]
         0   , 20;      % ddf, Degree-day factor for snowmelt [mm/oC/d]
         1   , 2000;    % Sb1, Maximum soil moisture storage [mm]
         0   , 1 ;      % tw, Groundwater leakage time [d-1]
         0   , 1 ;      % I_alpha, Intercepted fraction of Pr [-]
         1   , 365;     % I_s, Maximum Leaf Area Index timing [d]
         0   , 1 ;      % tu, Slow flow routing response time [d-1]
         0.05, 0.95;    % se, Root zone storage capacity as fraction of Sb2 [-]
         1   , 2000;    % Sb2, Root zone storage capacity [mm]
         0   , 1];      % tc, Mean residence time [d-1]

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'tcrit, Snowfall & snowmelt temperature [oC]' ...
                           'ddf, Degree-day factor for snowmelt [mm/oC/d]' ...
                           'Sb1, Maximum soil moisture storage [mm]' ...
                           'tw, Groundwater leakage time [d-1]' ...
                           'I_alpha, Intercepted fraction of Pr [-]' ...
                           'I_s, Maximum Leaf Area Index timing [d]' ...
                           'tu, Slow flow routing response time [d-1]' ...
                           'se, Root zone storage capacity as fraction of Sb2 [-]' ...
                           'Sb2, Root zone storage capacity [mm]' ...
                           'tc, Mean residence time [d-1]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)