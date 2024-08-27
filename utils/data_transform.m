function [RijMat,Ind,W,W_gt,AdjMat,G, R_gt] = data_transform(model)
AdjMat = model.AdjMat;
n = size(AdjMat,1);
Ind = model.Ind;
RijMat = model.RijMat;
R_orig = model.R_orig;
noiseInd = model.noiseInd;
corrInd = model.corrInd;
Ind_i = Ind(:,1);
Ind_j = Ind(:,2);
R_gt = zeros(3*n,3);
W = zeros(3*n,3*n);
G = zeros(n,n);
for i=1:n
    R_gt(3*i-2:3*i,:) = R_orig(:,:,i);
end
W_gt = R_gt*R_gt';
for l=1:length(Ind_i)
    i = Ind_i(l);j = Ind_j(l);
    W(3*i-2:3*i,3*j-2:3*j) = RijMat(:,:,l);
    W(3*j-2:3*j,3*i-2:3*i) = RijMat(:,:,l)';
end
for l=noiseInd
    i = Ind_i(l);j = Ind_j(l);
    G(i,j) = 1;G(j,i) = 1;
end
for l=corrInd
    i = Ind_i(l);j = Ind_j(l);
    G(i,j) = 0;G(j,i) = 0;
end