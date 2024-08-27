%% Author: Shaohan Li
%% © Regents of the University of Minnesota. All rights reserved
%%------------------------------------------------
%% Rotation Synchronization by IRLS using LongSync minimum spanning tree for initialization
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% RijMat: 3 by 3 by edge_num tensor that stores the given relative rotations corresponding to Ind
%% n: number of cameras
%% params.LongSync_option: 'only4' if using only 4-cycles, 'only5' if using only 5-cycles, '34' if using 3,4-cycles, '345' if using 3,4,5-cycles, '3' if using 3-cycles only
%% 
%% Output:
%% R_est: R_est: Estimated rotations (3 by 3 by n)
%% Weights_LongSync: Estimated edge weight matrix (n by n)

function [R_est,Weights_LongSync] = ...
    LongSync_IRLS(Ind,RijMat,n,params)
%% Build adjacency and stacked rotation matrices
[RijMat,Ind,W,AdjMat] = data_transform_real(Ind,RijMat,n);
%% Compute optimal weights
d = 3;
option = params.LongSync_option;
[Weights] = LongSync(Ind,n,W,d,option);
Weights_LongSync = Weights;

RijMat1 = RijMat;
Ind1 = Ind;

% Initialize with Minimum Spanning Tree
G_LongSync = graph((2-Weights_LongSync).*AdjMat,'upper');
C = conncomp(G_LongSync);
NodeLog = true(1,length(C));
G_LongSync = subgraph(G_LongSync,NodeLog);
Tree = minspantree(G_LongSync);
AdjTree = adjacency(Tree);
rootnodes = 1;
NodeMap = cumsum(NodeLog);
n = sum(NodeLog,'all');
IndLog2 = true(size(Ind1,1),1);
for i=1:size(Ind1,1)
    if NodeLog(Ind1(i,1)) && NodeLog(Ind1(i,2))
        Ind1(i,:) = NodeMap(Ind1(i,:));
    else
        IndLog2(i)=false;
    end
end
RijMat1 = RijMat1(:,:,IndLog2);
Ind1 = Ind1(IndLog2,:);
AdjMat = AdjMat(NodeLog,NodeLog);

if isempty(Ind1)
    Ind_i = zeros(0,1);
    Ind_j = zeros(0,1);
else
    Ind_i = Ind1(:,1);
    Ind_j = Ind1(:,2);
end
added=zeros(1,n);
R0 = zeros(3*n,3);
R0(3*rootnodes-2:3*rootnodes,:)=eye(3);
added(rootnodes)=1;
newroots = [];

IndMat = zeros(n,n);
for ll=1:length(Ind_i)
    i = Ind_i(ll);j = Ind_j(ll);
    IndMat(i,j) = ll;
    IndMat(j,i) = -ll;
end

iter=1;
while sum(added)<n && iter < 100
    for node_root = rootnodes
        leaves = find((AdjTree(node_root,:).*(1-added))==1);
        newroots = [newroots, leaves];
        for node_leaf=leaves
            edge_leaf = IndMat(node_leaf,node_root);
            if edge_leaf>0
                R0(3*node_leaf-2:3*node_leaf,:)=...
                    RijMat1(:,:,abs(edge_leaf))*R0(3*node_root-2:3*node_root,:);
            else
                R0(3*node_leaf-2:3*node_leaf,:)=...
                    (RijMat1(:,:,abs(edge_leaf)))'*R0(3*node_root-2:3*node_root,:);
            end
            added(node_leaf)=1;
        end
    end
    rootnodes = newroots;
    iter = iter+1;
end

R_est = zeros(3,3,n);
for i=1:n
    if trace(R0(3*i-2:3*i,:))<-1
        R0(3*i-2:3*i,:) = proj_SO3(R0(3*i-2:3*i,:));
    end
    R_est(:,:,i) = R0(3*i-2:3*i,:);
end

% Start IRLS

tmax = 300;
R_init_LS = R_est;
Ind_T = Ind1';
Amatrix = Build_Amatrix(Ind_T,n); % A matrix for least squares solver (by AVISHEK CHATTERJEE)      
RR = permute(RijMat1, [2,1,3]); % relative rotations -- take transpose as the original LAA code estimates R' in our setting
Q = R2Q(R_init_LS); % Transfer to quoternion representation (by AVISHEK CHATTERJEE)  
QQ = R2Q(RR); % Transfer to quoternion representation (by AVISHEK CHATTERJEE) 
Weights1_LS = Weights_LongSync;

tau_mpls = [1,1,1,1];
T_tau_mpls = length(tau_mpls);
if T_tau_mpls<tmax
    tau_mpls = [tau_mpls,tau_mpls(end)*(ones(1,tmax-T_tau_mpls))];
end
stop_threshold = 1e-3;
score = inf; Iteration = 1;
maxIters = 100;
weight_max = 1e4;
disp('Rotation Initialized!')
disp('Start IRLS reweighting ...')
% start IRLS reweighting iterations
while((score>stop_threshold)&&(Iteration<maxIters))
     % one iteration of Weighted Lie-Algebraic Averaging (by AVISHEK CHATTERJEE)
    if Iteration == 1
        Weights1_LS = Weights1_LS(sub2ind(size(Weights1_LS),Ind_i,Ind_j));
    end
    [Q,W,B,score] = Weighted_LAA(Ind_T,Q,QQ,Amatrix,Weights1_LS,n); 
    E=(Amatrix*W(2:end,2:4)-B); 
    ResVec = sqrt(sum(E.^2,2))/pi; % normalized residual rij for all edges
    RHVec = ResVec; % convex combination of rij and hij     
    Weights1_LS = (1./(RHVec.^0.75)); % compute edge weights
    % additional truncation for edge weights
    Weights1_LS(Weights1_LS>weight_max)= weight_max;
    % report the change of estimated rotations (stop when the change is small) 
    fprintf('Iter %d: ||\x394R||= %f\n', Iteration, score); 
    Iteration = Iteration+1;        
end
% transform from quaternion and return the estimated rotations
R_est=zeros(3,3,n);
for i=1:n
    if i<=size(Q,1)
        R_est(:,:,i)=q2R(Q(i,:));
    else
        R_est(:,:,i)=eye(3);
    end
    if sum(isnan(R_est(:,:,i)),'all')>0
        R_est(:,:,i) = eye(3);
    end 
end 
disp('DONE!')


end


function[Q] = genrandSO3()
    Q=randn(3);
    [U, ~, V]= svd(Q);
    S0 = diag([1,1,det(U*V')]);  
    Q=U*S0*V';
end

function [y] = g_tilde1(x,beta)
y = exp(-beta*sqrt(max(1-x,0)));
y = max(y,0);
y = min(y,1);
end

function [Y] = g_tilde(X,beta)
[m,n] = size(X);
Y = zeros(m,n);
[p,q,v] = find(X);
for t=1:length(p)
    i = p(t);
    j = q(t);
    Y(i,j) = g_tilde1(v(t),beta);
end
end

function [Q] = proj_SO3(Q)
    [U, ~, V]= svds(Q);
    S0 = diag([1,1,det(U*V')]);  
    Q=U*S0*V';
end

function out = sum_blocks(A,block_nrows, block_ncols)
out = squeeze(sum(reshape(sum(reshape(A,block_nrows,[])),...
                    size(A,1)/block_nrows,block_ncols,[]),2));
end

function [R0] = Rmat_to_Rmat3D(R)
[m,n] = size(R);
if mod(m,n)==0
    R0 = zeros(n,n,m/n);
    for i=1:m/n
        R0(:,:,i) = R(n*i-n+1:n*i,:);
    end
else
    error('Rmat_to_Rmat3D error: require input matrix size [m,n] with mod(m,n)==0.');
end
end

function [R0] = Rmat3D_to_Rmat(R)
[m,n,p] = size(R);
if m==n
    R0 = zeros(m*p,m);
    for i=1:p
        R0(m*i-m+1:m*i,:) = R(:,:,i);
    end
else
    error('Rmat3D_to_Rmat error: require input matrix size [m,n] with mod(n,m)==0.');
end
end