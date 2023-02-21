function [func] = interflow_12(~)
% interflow_12 
% source interflow_3 with customization -> fiels capacity (FC) is lower
% bounrdy 
%
% Anonymous function
% ------------------
% Description:  Non-linear interflow(variant) when current storage is over 
%               a threshhold (FC) and zero otherwise
% Constraints:  f <= S-FC 
%               S >= 0  - this avoids numerical issues with complex numbers
%               p2*Smax = FC
% @(Inputs):    p1   - time delay [d-1]
%               p2   - field capacity coefficient[-]
%               p3   - exponential scaling parameter [-]
%               S    - current storage [mm]
%               Smax - maximum storage [mm]
%               dt   - time step size [d]
%CB


func = @(p1,p2,p3,S,Smax,dt) ((S >(p2*Smax)).*(min(p1*max((S-(p2*Smax)),0).^(p3),max(S/dt,0))));

end

