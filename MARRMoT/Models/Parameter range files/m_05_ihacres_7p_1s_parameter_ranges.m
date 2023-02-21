function [ theta ] = m_05_ihacres_7p_1s_parameter_ranges( )
%m_05_ihacres_7p_1s_parameter_ranges Provides parameter ranges for calibration
%   of the IHACRES model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Croke, B. F. W., & Jakeman, A. J. (2004). A catchment moisture deficit 
% module for the IHACRES rainfall-runoff model. Environmental Modelling and
% Software, 19(1), 1�5. http://doi.org/10.1016/j.envsoft.2003.09.001
%
% Littlewood, I. G., Down, K., Parker, J. R., & Post, D. A. (1997). IHACRES
% v1.0 User Guide.
% 
% Ye, W., Bates, B. C., Viney, N. R., & Sivapalan, M. (1997). Performance 
% of conceptual rainfall-runoff models in low-yielding ephemeral 
% catchments. Water Resources Research, 33(1), 153�166. 
% http://doi.org/doi:10.1029/96WR02840

theta = [1   , 2000;    % lp, Wilting point [mm]
         1   , 2000;    % d, Threshold for flow generation [mm]
         0   , 10;      % p, Flow response non-linearity [-]
         0   , 1;       % alpha, Fast/slow flow division [-]
         1   , 700;     % tau_q, Fast flow routing delay [d]
         1   , 700;     % tau_s, Slow flow routing delay [d]
         0   , 119];    % tau_d, flow delay [d]

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'lp, Wilting point [mm]' ...
                           'd, Threshold for flow generation [mm]' ...
                           'p, Flow response non-linearity [-]' ...
                           'alpha, Fast/slow flow division [-]' ...
                           'tau_q, Fast flow routing delay [d]' ...
                           'tau_s, Slow flow routing delay [d]' ...
                           'tau_d, flow delay [d]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)     