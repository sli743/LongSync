%% Author: Shaohan Li
%% © Regents of the University of Minnesota. All rights reserved
%%------------------------------------------------
%% LongSync
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% n: number of cameras
%% R: d*n by d*n matrix that stores the given relative SO(d) measurements
%% For rotation synchronization, W(3*i-2:3*i, 3*j-2:3*j) stores the rotation between camera i and camera j
%% d: dimension of SO(d) matrix to be synchronized, d = 3 for Rotation Synchronization
%% option: string that indicates what cycle lengths to be used
%%      'only4' if using only 4-cycles
%%      'only5' if using only 5-cycles
%%      '34' if using 3,4-cycles
%%      '345' if using 3,4,5-cycles
%%      '3' if using 3-cycles only
%% Output:
%% Weights: Estimated edge weights (n by n)
%% SVec: Estimated edge corruption levels (1 by edge_num)

function [Weights,SVec] = LongSync(Ind,n,R,d,option)
if isempty(Ind)
    Ind_i = zeros(0,1);
    Ind_j = zeros(0,1);
else
    Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
end

beta = 1;
beta_rate = 2;
beta_max = 20;
max_iter = 10;

% Initialize weights 
AdjMat = sparse(Ind_i,Ind_j,ones(1,length(Ind_i)),n,n);
AdjMat = AdjMat + AdjMat';
Weights = tril(ones(n,n) .* AdjMat,-1);
Weights = sparse(Weights + Weights');
eps = 1e-8;
llist = [2,3,4];
ll = length(llist);
% w = [1,cind];% w is the indicator for using 3, 4 or 5-cycle information
switch true
    case strcmp(option, 'only4')
        w = [0,1,0];
    case strcmp(option, 'only5')
        w = [0,0,1];
    case strcmp(option, '3')
        w = [1,0,0];
    case strcmp(option, '34')
        w = [1,1,0];
    case strcmp(option, '345')
        w = [1,1,1];
end
AdjMat1 = AdjMat;
thres_w = 0;
Adj = cell(1,ll);
for i=1:ll
    Adj{i} = zeros(n,n);
end
A = AdjMat1;
a3 = A^3.*A;
a2 = A^2.*A>thres_w;
degVec = sum(A);
degMat = (degVec+degVec'-1).*A;
Adj{1} = a2.*A;
Adj{2} = (a3-degMat>thres_w).*(A-Adj{1});
AdjMat2 = Adj{1}+w(2)*Adj{2};
if w(3)==1
    a4 = A^4-(degVec+degVec'-3*A).*(A^2)-diag(diag(A^3))*A-A*diag(diag(A^3))-A*diag(degVec-2)*A;
    a4 = a4>thres_w;
    Adj{3} = a4.*(A-Adj{1}-Adj{2});
    AdjMat2 = AdjMat2 + Adj{3};
end
Adjbad = AdjMat1 - AdjMat2;
Weights = Weights.*AdjMat2;
if strcmp(option, 'only4')
    Adj{2} = AdjMat2;
    Adj{1} = sparse(n,n);
    Adj{3} = sparse(n,n);
    max_iter = 100;
end
if strcmp(option, 'only5')
    Adj{3} = AdjMat2;
    Adj{1} = sparse(n,n);
    Adj{2} = sparse(n,n);
    max_iter = 100;
end
for iter = 1:max_iter
    P1 = zeros(d*n,d*n);
    P2 = zeros(d*n,d*n);
    c = 3;
    X2 = gc(Weights,R,d,c);
    Y2 = fc(Weights,c);
    P1 = P1 + X2.*kron(Adj{1},ones(d,d));
    P2 = P2 + kron(Y2,ones(d,d)).*kron(Adj{1},ones(d,d));
    if w(2) == 1
        c = 4;
        X3 = gc(Weights,R,d,c);
        Y3 = fc(Weights,c);
        P1 = P1 + X3.*kron(Adj{2},ones(d,d));
        P2 = P2 + kron(Y3,ones(d,d)).*kron(Adj{2},ones(d,d));
    end
    if w(3) == 1
        c = 5;
        X4 = gc(Weights, R,d,c);
        Y4 = fc(Weights,c);
        P1 = P1 + X4.*kron(Adj{3},ones(d,d));
        P2 = P2 + kron(Y4,ones(d,d)).*kron(Adj{3},ones(d,d));
    end
    P0 = spdd(P1,P2,-1);
    P = P0.*R;
    P = sconv2(P,sparse(ones(d,d)));
    P = P(d*(1:n),d*(1:n))/d;
    P = min(P,1);P = max(P,-1/d);
    Weights_0 = Weights;
    tic;
    Weights = g_tilde(P,beta).* AdjMat2;
    diff = norm(Weights_0-Weights,'F')^2/norm(Weights_0,'F')^2;
    if diff < eps
        break;
    end
    beta = min(beta*beta_rate,beta_max);
end
Weights(Adjbad>0) = exp(-beta);
ne = length(Ind_i);
SVec = zeros(1,ne);
for k=1:ne
    i = Ind_i(k);
    j = Ind_j(k);
    SVec(k) = sqrt(1-P(i,j));
end
end

function [y] = g_tilde1(x,beta)
y = exp(-beta*sqrt(max(1-x,0)));
y = max(y,0);
y = min(y,1);
end

function [Y] = g_tilde(X,beta)
[m,n] = size(X);
[p,q,v] = find(X);
Y = zeros(m,n);
for k=1:length(v)
    i = p(k);j = q(k);
    Y(i,j) = g_tilde1(v(k),beta);
end
end

function [X3] = gc(Weights,W,d,c)
if c==3
    X3 = (kron(Weights,ones(d,d)) .* W)^2;
end
if c==4
    Y3w = Weights.*diag(Weights^2)+Weights.*diag(Weights^2)'-Weights.^3;
    Y3 = kron(Y3w,ones(d,d));
    X3 = (Y3 .* W);
    X3 = (kron(Weights,ones(d,d)) .* W)^3 - X3;
end
if c==5
    n = size(Weights,1);
    M = kron(eye(n), ones(d,d));
    P = kron(Weights,ones(d,d)).*W;
    W2k = kron(Weights.^2,ones(d,d));
    blkT = @(P) blockproc(full(P),[3,3],@(Q) Q.data');
    blkmul = @(P,Q) cell2mat(cellfun(@(A,B) A*B, mat2cell(P,d*ones(1,n),...
        d*ones(1,n)),mat2cell(Q,d*ones(1,n),d*ones(1,n)),'UniformOutput',false));
    P2 = P^2;
    P3 = P^3;
    X3 = P^4 - (P3.*M)*P - (P2.*M)*P2 - P2*(P2.*M) - ...
        P*(P3.*M) - P*(P2.*M)*P + sparse(blkmul(blkmul(P, blkT(P2)),P)) +...
        W2k.*P2 + P2.*W2k + P*(W2k.*P) + (W2k.*P)*P;
end
end

function [Y3] = fc(Weights,c)
if c==3
    Y3 = Weights^2;
end
if c==4
    Y3w = Weights.*diag(Weights^2)+Weights.*diag(Weights^2)'-Weights.^3;
    Y3 = Weights^3 - Y3w;
end
if c==5
    Y3 = Weights^4;
    W3 = Weights^3;
    W2 = Weights^2;
    W1 = Weights;
    Y3 = Y3 - W1.*diag(W3) - W2.*diag(W2) - W2.*diag(W2)' - W1.*diag(W3)'...
        - W1*(W1.*diag(W2)) + 3*(W1.*W1).*W2 + W1*(W1.*W1.*W1)...
        + (W1.*W1.*W1)*W1;
end
end
