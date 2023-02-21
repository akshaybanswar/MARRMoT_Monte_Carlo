function [func] = infiltration_8(~)
%infiltration_8
%
% Anonymous function
% ------------------
% Description:  Infiltration into storage is equal the inflow when current 
%               storage is under the maximum storage, 
%               and zero when storage reaches maximum capacity 
% Constraints:  f <= fin
% @(Inputs):    
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               fin  - size of incoming flux [mm/d]
%CB

func = @(S,Smax,fin) ((S < Smax) .* fin);


end