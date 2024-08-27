addpath ./Models
median_error_LS_all = zeros(1,nr);
median_error_LS5_all = zeros(1,nr);
median_error_MPLS_all = zeros(1,nr);
median_error_IRLS_all = zeros(1,nr);
% median_error_LS_dist_all = zeros(1,nr);
% median_error_LS_dist53_all = zeros(1,nr);
mean_error_LS_all = zeros(1,nr);
mean_error_LS5_all = zeros(1,nr);
mean_error_MPLS_all = zeros(1,nr);
mean_error_IRLS_all = zeros(1,nr);
% mean_error_LS_dist_all = zeros(1,nr);
% mean_error_LS_dist53_all = zeros(1,nr);
for r = 1:nr
rng(r);
% q = 0.88;%for adv set q < 0.5, for unif set q < 1
% p0 = 1;
p1 = 0;
p2 = 1;
% q1 = 0.8;
% q2 = 0.8;
q1 = q; 
q2 = q;
n = 200;
% nb = [n/4,n/4,n/4,n/4];
% idx = [ones(1,n/4),2*ones(1,n/4),3*ones(1,n/4),4*ones(1,n/4)];
nb = [n/2,n/2];
idx = [ones(1,n/2),2*ones(1,n/2)];
d = 3;
sigma = 0;

if strcmp(option,'Unif')
    [RijMat,Ind,W,W_gt,AdjMat,G, R_gt] = Unif_corr_data_gen_SO3(q,p0,n);
end
if strcmp(option,'Adv')
    [RijMat,Ind,W,W_gt,AdjMat,G, R_gt] = Adv_corr_data_gen_SO3(q,p0,n);
end
if strcmp(option,'Unif_noise')
    [model_out]=Uniform_Topology(n,p0,q,sigma,'uniform');
    [RijMat,Ind,W,W_gt,AdjMat,G, R_gt] = data_transform(model_out);
end
if strcmp(option,'Bad_inlier_graph')
    [model_out]=Bad_Inlier_Graph_Topology(n,p1,p2,q,sigma,'Bad_inlier_graph');
    [RijMat,Ind,W,W_gt,AdjMat,G, R_gt] = data_transform(model_out);
end
if strcmp(option,'Sto_block')
    [model_out]=Sto_block(nb,p1,p2,q1,q2,sigma,'uniform');
    [RijMat,Ind,W,W_gt,AdjMat,G, R_gt] = data_transform(model_out);
end
G0 = max(G,0);
G0 = diag(sum(G0,2))\G0;

% 4-cycles only
R_gt3D = Rmat_to_Rmat3D(R_gt);
params.reweighting = 32;
params.cycle_indicator = [1,0];
t0 = cputime;
params.LongSync_option = 'only4';
[R_LS3D] = ...
    LongSync_IRLS(Ind,RijMat,n,params, 0, 0, R_gt3D, params.cycle_indicator);
t_LS = cputime-t0;
% R_LS3D = Rmat_to_Rmat3D(R0);
fprintf('Mean/Median Error By Rotation Alignment:\n')
if strcmp(option_align, 'L2')
    [R_out_LS3D, R_align_LS3D, mean_error_LS, median_error_LS]...
        = Rotation_Alignment(R_LS3D, R_gt3D);
else
    error_LS= error_R(R_LS3D, R_gt3D);
    mean_error_LS = mean(error_LS);
    median_error_LS = median(error_LS);
end

% 5-cycles only
R_gt3D = Rmat_to_Rmat3D(R_gt);
params.reweighting = 32;
params.cycle_indicator = [1,1];
t0 = cputime;
params.LongSync_option = 'only5';
[R_LS3D5] = ...
    LongSync_IRLS(Ind,RijMat,n,params, 0, 0, R_gt3D, params.cycle_indicator);
t_LS5 = cputime-t0;
% R_LS3D = Rmat_to_Rmat3D(R0);
fprintf('Mean/Median Error By Rotation Alignment:\n')
% [R_out_LS3D5, R_align_LS3D5, mean_error_LS5, median_error_LS5]...
%     = Rotation_Alignment(R_LS3D5, R_gt3D);
if strcmp(option_align, 'L2')
    [R_out_LS3D5, R_align_LS3D5, mean_error_LS5, median_error_LS5]...
        = Rotation_Alignment(R_LS3D5, R_gt3D);
else
    error_LS5= error_R(R_LS3D5, R_gt3D);
    mean_error_LS5 = mean(error_LS5);
    median_error_LS5 = median(error_LS5);
end

t0 = cputime;
CEMP_parameters.max_iter = 6;
CEMP_parameters.reweighting = 2.^((1:6)-1);
CEMP_parameters.nsample = 50;

% set MPLS parameters
MPLS_parameters.stop_threshold = 1e-3;
MPLS_parameters.max_iter = 100;
% MPLS_parameters.thresholding = [0.95,0.9,0.85,0.8];
MPLS_parameters.thresholding = [1,1,1,1];
MPLS_parameters.reweighting = CEMP_parameters.reweighting(end);
% MPLS_parameters.cycle_info_ratio = 1./((1:MPLS_parameters.max_iter)+1);
MPLS_parameters.cycle_info_ratio = 1./((1:MPLS_parameters.max_iter)+1);
[R_MPLS] = MPLS(Ind,RijMat,CEMP_parameters, MPLS_parameters);
t_MPLS = cputime-t0;
R_MPLS3D = R_MPLS;



if strcmp(option_align, 'L2')
    [R_out_MPLS3D, R_align_MPLS3D, mean_error_MPLS, median_error_MPLS]...
        = Rotation_Alignment(R_MPLS3D, R_gt3D);
else
    error_MPLS= error_R(R_MPLS3D, R_gt3D);
    mean_error_MPLS = mean(error_MPLS);
    median_error_MPLS = median(error_MPLS);
end


t0 = cputime;
CEMP_parameters.max_iter = 6;
CEMP_parameters.reweighting = 2.^((1:6)-1);
CEMP_parameters.nsample = 50;

% set MPLS parameters
MPLS_parameters.stop_threshold = 1e-3;
MPLS_parameters.max_iter = 100;
% MPLS_parameters.thresholding = [0.95,0.9,0.85,0.8];
MPLS_parameters.thresholding = [1,1,1,1];
MPLS_parameters.reweighting = CEMP_parameters.reweighting(end);
% MPLS_parameters.cycle_info_ratio = 1./((1:MPLS_parameters.max_iter)+1);
MPLS_parameters.cycle_info_ratio = 0;
[R_IRLS] = IRLS(Ind,RijMat,CEMP_parameters, MPLS_parameters);
t_IRLS = cputime-t0;
R_IRLS3D = R_IRLS;



if strcmp(option_align, 'L2')
    [R_out_IRLS3D, R_align_IRLS3D, mean_error_IRLS, median_error_IRLS]...
        = Rotation_Alignment(R_IRLS3D, R_gt3D);
else
    error_IRLS= error_R(R_IRLS3D, R_gt3D);
    mean_error_IRLS = mean(error_IRLS);
    median_error_IRLS = median(error_IRLS);
end

mean_error_LS
mean_error_LS5
mean_error_MPLS
median_error_LS
median_error_LS5
median_error_MPLS
% mean_error_LS_dist
% median_error_LS_dist
% mean_error_LS_dist53
% median_error_LS_dist53


median_error_LS_all(r) = median_error_LS;
median_error_LS5_all(r) = median_error_LS5;
median_error_MPLS_all(r) = median_error_MPLS;
median_error_IRLS_all(r) = median_error_IRLS;
% median_error_LS_dist_all(r) = median_error_LS_dist;
% median_error_LS_dist53_all(r) = median_error_LS_dist53;
mean_error_LS_all(r) = mean_error_LS;
mean_error_LS5_all(r) = mean_error_LS5;
mean_error_MPLS_all(r) = mean_error_MPLS;
mean_error_IRLS_all(r) = mean_error_IRLS;
% mean_error_LS_dist_all(r) = mean_error_LS_dist;
% mean_error_LS_dist53_all(r) = mean_error_LS_dist53;
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