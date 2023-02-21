function [ theta ] = m_40_smar_8p_6s_parameter_ranges( )
%m_40_smar_8p_6s_parameter_ranges Provides parameter ranges for calibration
%   of the SMAR model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%
% Model reference
% O�Connell, P. E., Nash, J. E., & Farrell, J. P. (1970). River flow 
% forecasting through conceptual models part II - the Brosna catchment at 
% Ferbane. Journal of Hydrology, 10, 317�329.
%
% Tan, B. Q., & O�Connor, K. M. (1996). Application of an empirical 
% infiltration equation in the SMAR conceptual model. Journal of Hydrology,
% 185(1-4), 275�295. http://doi.org/10.1016/0022-1694(95)02993-1

theta = [   0,   1;     % h, Maximum fraction of direct runoff [-] 
            0, 200;     % y, Infiltration rate [mm/d] 
            1,2000;     % smax, Maximum soil moisture storage [mm]
            0,   1;     % c, Evaporation reduction coefficient [-]
            0,   1;     % g, Groundwater recharge coefficient [-]
            0,   1;     % kg, Groundwater time parameter [d-1]
            1,  10;     % n, Number of Nash cascade reservoirs [-]
            1, 120];    % n*k, Routing delay [d]
      
% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'h, Maximum fraction of direct runoff [-] ' ...
                           'y, Infiltration rate [mm/d] ' ...
                           'smax, Maximum soil moisture storage [mm]' ...
                           'c, Evaporation reduction coefficient [-]' ...
                           'g, Groundwater recharge coefficient [-]' ...
                           'kg, Groundwater time parameter [d-1]' ...
                           'n, Number of Nash cascade reservoirs [-]' ...
                           'n*k, Routing delay [d]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)