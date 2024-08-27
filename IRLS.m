%% Author: Shaohan Li
%% © Regents of the University of Minnesota. All rights reserved
%%------------------------------------------------
%% IRLS for Rotation Synchronization
%%------------------------------------------------
%% Input Parameters: 
%% Ind: edge_num by 2 "edge indices matrix". Each row is the index of an edge (i,j). that is sorted as (1,2), (1,3), (1,4),... (2,3), (2,4),.... 
%% edge_num is the number of edges.
%% RijMat: 3 by 3 by edge_num tensor that stores the given relative rotations corresponding to Ind

%% CEMP_parameters.max_iter: total # of iterations for CEMP
%% CEMP_parameters.reweighting: the sequence of reweighting parameter beta_t for CEMP
%% CEMP_parameters.nsample: # of cycles sampled per edge (also used in MPLS stage)

%% MPLS_parameters.stop_threshold: stopping criterion
%% MPLS_parameters.max_iter: the maximal number of iterations of MPLS
%% MPLS_parameters.reweighting: the sequence of reweighting parameter beta_t for MPLS
%% MPLS_parameters.thresholding: the sequence of thresholding parameters tau_t
%% MPLS_parameters.cycle_info_ratio: the coefficient alpha_t of cycle-consistency information.



%% Output:
%% R_est: Estimated rotations (3x3xn)
%% R_init: Initialized rotations by CEMP+MST (3x3xn)

%% Reference
%% [1] Yunpeng Shi and Gilad Lerman. "Message Passing Least Squares Framework and its Application to Rotation Synchronization" ICML 2020.


function[R_est, R_init, Weights, SVec] = IRLS(Ind,RijMat,CEMP_parameters, MPLS_parameters, R_init)

    %CEMP parameters
    T=CEMP_parameters.max_iter; 
    beta_cemp=CEMP_parameters.reweighting;
    nsample = CEMP_parameters.nsample;
    T_beta = length(beta_cemp);
    if T_beta<T
        % if the reweighting parameter vector is short, then the rest of
        % missing elements are set to constant
        beta_cemp = [beta_cemp,beta_cemp(end)*(ones(1,T-T_beta))]; 
    end
    
    % MPLS paramters
    stop_threshold=MPLS_parameters.stop_threshold;
    maxIters = MPLS_parameters.max_iter;    
    % building the graph   
    Ind_i = Ind(:,1);
    Ind_j = Ind(:,2);
    if isfield(MPLS_parameters,'n')
        n = MPLS_parameters.n;
    else
        n=max(Ind,[],'all');
    end
    m=size(Ind_i,1);
        
    % start IRLS procedure
    
    % transform the data format for the following Lie-Alegbraic
    % Averaging (LAA) solver that was written by AVISHEK CHATTERJEE
    RR = permute(RijMat, [2,1,3]); % relative rotations -- take transpose as the original LAA code estimates R' in our setting
    Ind_T = Ind'; % indices matrix
    % Formation of A matrix.
    Amatrix = Build_Amatrix(Ind_T,n); % A matrix for least squares solver (by AVISHEK CHATTERJEE)      
    Q = R2Q(R_init); % Transfer to quoternion representation (by AVISHEK CHATTERJEE)  
    QQ = R2Q(RR); % Transfer to quoternion representation (by AVISHEK CHATTERJEE)  
    score=inf;    Iteration=1;
 
    % initialization edge weights using (sij) estimated by CEMP
    Weights = ones(m,1);
    weight_max = 1e4; % weights cannot be too large nor too small (for numerical stability and graph connectivity)
    Weights(Weights>weight_max)= weight_max;
   
    disp('Rotation Initialized!')
    disp('Start MPLS reweighting ...')
    % start MPLS reweighting iterations
    while((score>stop_threshold)&&(Iteration<maxIters))
         % one iteration of Weighted Lie-Algebraic Averaging (by AVISHEK CHATTERJEE)                
        [Q,W,B,score] = Weighted_LAA(Ind_T,Q,QQ,Amatrix,Weights,n); 
        E=(Amatrix*W(2:end,2:4)-B); 
        ResVec = sqrt(sum(E.^2,2))/pi; % normalized residual rij for all edges
        RHVec = ResVec; %  Only use residual for reweighting
        Weights = (1./(RHVec.^0.75)); % compute edge weights
        % additional truncation for edge weights
        Weights(Weights>weight_max)= weight_max;
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


function [Q] = proj_SO3(Q)
    [U, ~, V]= svds(Q);
    S0 = diag([1,1,det(U*V')]);  
    Q=U*S0*V';
end