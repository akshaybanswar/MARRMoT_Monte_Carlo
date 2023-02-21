function [ theta ] = m_20_gsfb_8p_3s_parameter_ranges( )
%m_20_gsfb_8p_3s_parameter_ranges Provides parameter ranges for calibration
%   of the 3-store GSFB model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Nathan, R. J., & McMahon, T. A. (1990). SFB model part l . Validation of 
% fixed model parameters. In Civil Eng. Trans. (pp. 157�161).
% 
% Ye, W., Bates, B. C., Viney, N. R., & Sivapalan, M. (1997). Performance 
% of conceptual rainfall-runoff models in low-yielding ephemeral catchments.
% Water Resources Research, 33(1), 153�166. http://doi.org/doi:10.1029/96WR02840

theta = [   0, 1;           % c, Recharge time coeffcient [d-1]
            0.05, 0.95;     % ndc, Threshold fraction of Smax [-]
            1, 2000;        % smax, Maximum soil moisture storage [mm]
            0, 20;          % emax, Maximum evaporation flux [mm/d]
            0, 200;         % frate, Maximum infiltration rate [mm/d]
            0, 1;           % b, Fraction of subsurface flow that is baseflow [-]
            0, 1;           % dpf, Baseflow time coefficient [d-1]
            1, 300];        % sdrmax, Threshold before baseflow can occur [mm]

   
% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'c, Recharge time coeffcient [d-1]' ...
                           'ndc, Threshold fraction of Smax [-]' ...
                           'smax, Maximum soil moisture storage [mm]' ...
                           'emax, Maximum evaporation flux [mm/d]' ...
                           'frate, Maximum infiltration rate [mm/d]' ...
                           'b, Fraction of subsurface flow that is baseflow [-]' ...
                           'dpf, Baseflow time coefficient [d-1]' ...
                           'sdrmax, Threshold before baseflow can occur [mm]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)