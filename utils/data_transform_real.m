function [RijMat,Ind,W,AdjMat] = data_transform_real(Ind,RijMat,n)
if isempty(Ind)
    AdjMat = sparse(n,n);
    W = zeros(3*n,3*n);
else
    AdjMat = sparse(Ind(:,1),Ind(:,2),ones(size(Ind,1),1),n,n);
    AdjMat = AdjMat+AdjMat';AdjMat = min(AdjMat,1);
    W = zeros(3*n,3*n);
    for l=1:size(Ind,1)
        i = Ind(l,1);j = Ind(l,2);
        W(3*i-2:3*i,3*j-2:3*j) = RijMat(:,:,l);
        W(3*j-2:3*j,3*i-2:3*i) = RijMat(:,:,l)';
    end
end
end