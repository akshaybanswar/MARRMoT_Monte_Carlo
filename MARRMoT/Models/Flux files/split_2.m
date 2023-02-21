function [func] = split_2(~)
%split_2 Creates function for flow splitting, counterpart to split_1
%
%
% Anonymous function
% ------------------
% Description:  Split flow (returns flux [mm/d])
% Constraints:  -
% @(Inputs):    p1   - fraction of flux to be diverted [-]
%               In   - incoming flux [mm/d]
%
%CB

func = @(p1,In) (1-p1).*In;

end