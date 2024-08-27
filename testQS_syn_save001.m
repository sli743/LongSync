addpath ./Models
median_error_LS_all = zeros(1,nr);
median_error_LS5_all = zeros(1,nr);
median_error_MPLS_all = zeros(1,nr);
median_error_IRLS_all = zeros(1,nr);
mean_error_LS_all = zeros(1,nr);
mean_error_LS5_all = zeros(1,nr);
mean_error_MPLS_all = zeros(1,nr);
mean_error_IRLS_all = zeros(1,nr);
for r = 1:nr
rng(r);
p1 = 0;
p2 = 1;
q1 = q; 
q2 = q;
n = 200;
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
t0 = cputime;
params.LongSync_option = 'only4';
[R_LS3D] = ...
    LongSync_IRLS(Ind,RijMat,n,params);
t_LS = cputime-t0;
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
t0 = cputime;
params.LongSync_option = 'only5';
[R_LS3D5] = ...
    LongSync_IRLS(Ind,RijMat,n,params);
t_LS5 = cputime-t0;
fprintf('Mean/Median Error By Rotation Alignment:\n')
if strcmp(option_align, 'L2')
    [R_out_LS3D5, R_align_LS3D5, mean_error_LS5, median_error_LS5]...
        = Rotation_Alignment(R_LS3D5, R_gt3D);
else
    error_LS5= error_R(R_LS3D5, R_gt3D);
    mean_error_LS5 = mean(error_LS5);
    median_error_LS5 = median(error_LS5);
end

% MPLS
t0 = cputime;
CEMP_parameters.max_iter = 6;
CEMP_parameters.reweighting = 2.^((1:6)-1);
CEMP_parameters.nsample = 50;

% set MPLS parameters
MPLS_parameters.stop_threshold = 1e-3;
MPLS_parameters.max_iter = 100;
MPLS_parameters.thresholding = [1,1,1,1];
MPLS_parameters.reweighting = CEMP_parameters.reweighting(end);
MPLS_parameters.cycle_info_ratio = 1./((1:MPLS_parameters.max_iter)+1);
[R_MPLS, R_init] = MPLS(Ind,RijMat,CEMP_parameters, MPLS_parameters);
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

% IRLS
t0 = cputime;
CEMP_parameters.max_iter = 6;
CEMP_parameters.reweighting = 2.^((1:6)-1);
CEMP_parameters.nsample = 50;

% set MPLS parameters
MPLS_parameters.stop_threshold = 1e-3;
MPLS_parameters.max_iter = 100;
MPLS_parameters.thresholding = [1,1,1,1];
MPLS_parameters.reweighting = CEMP_parameters.reweighting(end);
MPLS_parameters.cycle_info_ratio = 0;
[R_IRLS] = IRLS(Ind,RijMat,CEMP_parameters, MPLS_parameters, R_init);
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
mean_error_IRLS
median_error_LS
median_error_LS5
median_error_MPLS
median_error_IRLS


median_error_LS_all(r) = median_error_LS;
median_error_LS5_all(r) = median_error_LS5;
median_error_MPLS_all(r) = median_error_MPLS;
median_error_IRLS_all(r) = median_error_IRLS;
mean_error_LS_all(r) = mean_error_LS;
mean_error_LS5_all(r) = mean_error_LS5;
mean_error_MPLS_all(r) = mean_error_MPLS;
mean_error_IRLS_all(r) = mean_error_IRLS;
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
