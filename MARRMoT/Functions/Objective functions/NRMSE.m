function r = NRMSE(Obs,Sim)
% INPUT 
% Obs M x N
% Sim M x N
% Output
% Result-struct

%% getting size and condition checking
[row_R,col_R,dim_R]=size(Obs);
[row_T,col_T,dim_T]=size(Sim);
if row_R~=row_T || col_R~=col_T || dim_R~=dim_T
    error('Input must have same dimensions')
end

%% Common function for matrix
% Mean for Matrix
meanmat=@(a)(mean(mean(a)));
% Sum for Matrix
summat=@(a)(sum(sum(a)));
% Min  for Matrix
minmat=@(a)(min(min(a)));
% Max  for Matrix
maxmat=@(a)(max(max(a)));

%% RMSE Root-mean-square deviation
rmse =abs( sqrt( meanmat((Sim-Obs).^2) ) );

%% Normalized RMSE Normalized Root-mean-square deviation
r =rmse/(maxmat(Obs)-minmat(Obs));
