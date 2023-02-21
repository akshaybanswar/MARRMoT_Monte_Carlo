function [func] = infiltration_7(~)
%infiltration_7
% source: infiltration_1 with customization
%
% Anonymous function
% ------------------
% Description:  Infiltration as exponentially declining based on relative storage
% Constraints:  f <= fin
% @(Inputs):    p1   - maximum infiltration rate [mm,/d]
%               p2   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               fin  - size of incoming flux [mm/d]
%CB

func = @(p1,p2,S,Smax,fin) ((S < Smax) .* (min(p1.*exp((-1*p2*S)./Smax),fin))).*(smoothThreshold_storage_logistic(S,Smax));


end