function r = pearson(Obs,Sim)
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
% Sum for Matrix
summat=@(a)(sum(sum(a)));

%% R Value
r =1-abs( summat((Sim-Obs).^2) / summat(Obs.^2) );
