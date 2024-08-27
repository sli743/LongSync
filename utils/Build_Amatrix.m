function[Amatrix] = Build_Amatrix(I,varargin)
m=size(I,2);
if length(varargin)>0
    N = varargin{1};
else
    N=max(max(I));
end
i=[[1:m];[1:m]];i=i(:);
j=I(:);
s=repmat([-1;1],[m,1]);
k=(j~=1);
Amatrix=sparse(i(k),j(k)-1,s(k),m,N-1);
end