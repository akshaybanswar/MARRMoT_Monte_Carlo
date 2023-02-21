function [ theta ] = m_36_modhydrolog_15p_5s_parameter_ranges( )
%m_36_modhydrolog_15p_5s_parameter_ranges Provides parameter ranges for calibration
%   of the MODHYDROLOG model.
%
% Copyright (C) 2018 W. Knoben
% This program is free software (GNU GPL v3) and distributed WITHOUT ANY
% WARRANTY. See <https://www.gnu.org/licenses/> for details.
%   
% Model reference
% Chiew, F. H. S. (1990). Estimating groundwater recharge using an
% integrated surface and groundwater model. University of Melbourne.
%
% Chiew, F., & McMahon, T. (1994). Application of the daily rainfall-runoff
% model MODHYDROLOG to 28 Australian catchments. Journal of Hydrology, 
% 153(1�4), 383�416. https://doi.org/10.1016/0022-1694(94)90200-3

theta = [0,     5;      % INSC, Maximum interception capacity, [mm]
         0,     600;    % COEFF, Maximum infiltration loss parameter, [mm]
         0,     15;     % SQ, Infiltration loss exponent, [-]
         1,     2000;   % SMSC, Maximum soil moisture capacity, [mm]
         0,     1;      % SUB, Proportionality constant, [-]
         0,     1;      % CRAK, Proportionality constant, [-]
         0,     20;     % EM, maximum plant-controled evap rate, [mm/d]
         0,     50;     % DSC, Maximum depression capacity, [mm]
         0,     1;      % ADS, Land fraction functioning as depression storage, [-]
         0.99,  1;      % MD, Depression storage parameter, [-]
         0,     0.5;    % VCOND, Leakage coefficient, [mm/d]
       -10,     10;     % DLEV, Datum around which groundwater fluctuates relative to river bed, [mm]
         0,     1;      % K1, Flow exchange parameter, [d-1] 
         0,     1;      % K2, Flow exchange parameter, [d-1] 
         0,     100];   % K3, Flow exchange parameter, [d-1] 

% Display an overview and warning      
txt = array2table(theta);
txt.Properties.VariableNames = {'min' 'max'};
txt.Properties.RowNames = {'INSC, Maximum interception capacity, [mm]' ...
                           'COEFF, Maximum infiltration loss parameter, [mm]' ...
                           'SQ, Infiltration loss exponent, [-]' ...
                           'SMSC, Maximum soil moisture capacity, [mm]' ...
                           'SUB, Proportionality constant, [-]' ...
                           'CRAK, Proportionality constant, [-]' ...
                           'EM, maximum plant-controled evap rate, [mm/d]' ...
                           'DSC, Maximum depression capacity, [mm]' ...
                           'ADS, Land fraction functioning as depression storage, [-]' ...
                           'MD, Depression storage parameter, [-]' ...
                           'VCOND, Leakage coefficient, [mm/d]' ...
                           'DLEV, Datum around which groundwater fluctuates relative to river bed, [mm]' ...
                           'K1, Flow exchange parameter, [d-1]' ...
                           'K2, Flow exchange parameter, [d-1]' ...
                           'K3, Flow exchange parameter, [d-1]'};
str = mfilename;                       
str(end-16:end) = [];

disp('---')
disp(['Overview of currently used parameter ranges for ',str,'.'])
disp('In the model parameter range files you can: (1) disable this warning, (2) adjust these ranges.')
disp('Please note that the MARRMoT default ranges are optional and not necessarily appropriate for your problem.')
disp(txt)     