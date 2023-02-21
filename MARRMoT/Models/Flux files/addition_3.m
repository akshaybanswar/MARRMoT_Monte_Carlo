function [func] = addition_3(~)
% addition_3
% 
% @(Inputs):    fin1   - incomming flux 1
%               fin2   - incomming flux 2
%               fin3   - incomming flux 3
%CB

func = @(fin1,fin2,fin3) (fin1 + fin2 + fin3);

end