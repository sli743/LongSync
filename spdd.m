function [A] = spdd(X,Y,s)
% Compute X./Y for nonzero entries only
% try extract nnz entries
[p,q] = size(X);
[i,j,v] = find(Y);
l = length(v);
v1 = zeros(l,1);
v2 = s*ones(l,1);
v1(v>0) = 1./v(v>0);
v2(v>0) = 0;
V1 = sparse(i,j,v1,p,q);
V2 = sparse(i,j,v2,p,q);
A = X.*V1-V2;
