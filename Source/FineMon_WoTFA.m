
%parameters:
% tensor size: I = time dimension, J = cycle dimension, K = metric dimension
% M = gpuArray(M);
I=size(M,1);
J=size(M,2);
K=size(M,3);
% window size after ST,  W_size = WT-T+1
W_size=101;
% \eta in residual(.)
% theta=2.8e-3;%MMS
% yita=4e-29;
% theta=2.4e-4;%SR
% yita=6e-28;
% theta=9e-4;%AQI
% yita=4e-29;  
% theta=2e-3;%CAIDA
% yita=8e-29;
theta=2.5e-3;%MAWI
yita=6e-29;
% theta=2e-4;%SMD
% yita=4e-29;
% \beta in delay-recovery strategy in Section 4.3
beta=3;
p=0.5;
% Whether to enable capturing large tubes by refining the frequency, the default is off because the refined frequency leads to an increase in the sample ratio
isRefine = 0;
% MMS:
% epsilon_delta = 2.7;
% epsilon_gamma = 0.2;
% CAIDA:
epsilon_delta = 2.8;
epsilon_gamma = 0.2;
% MAWI
% epsilon_delta = 2;
% epsilon_gamma = 0.7;

% perform the finemon for real dataset M
% [R,omega,h_incoms,estimators,rs,ms,com, pareto_oemga] = finemon(M,W_size,I,J,K,theta,yita,beta, isRefine, epsilon_delta, epsilon_gamma);
[R,omega,h_incoms,estimators,rs,ers, ms,com, pareto_oemga_large,pareto_oemga_low] = finemon(M,W_size,I,J,K,theta,yita,beta, isRefine, epsilon_delta, epsilon_gamma,p);            
% R = gether(R);
Rm=[];
for j=W_size+1:J
     Rm(size(Rm,1)+1:size(Rm,1)+I,:)=R(:,j,:);
end

% Performance Evaluation Criteria.
% inverse normalization
% N_max_min=MMS2_max_min;
% N_min=MMS2_min;
% N_max_min=KPIs2_max_min;
% N_min=KPIs2_min;
% N_max_min=Mone2_max_min;
% N_min=Mone2_min;
N_max_min=MAWI2_max_min;
N_min=MAWI2_min;
% N_max_min=CAIDA2_max_min;
% N_min=CAIDA2_min;
% N_max_min=SMD2_max_min;
% N_min=SMD2_min;
RM_orign=[];
for k=1:K
      RM_orign(:,k) = Rm(:,k)*N_max_min(k)+N_min(k);
end
% orig_data_M=MMS(8*40+1:(J-W_size+8)*40,:);%MMS
% orig_data_M=MMS(3*40+1:(J-W_size+3)*40,:);%MMS
% orig_data_M=KPIs(15*50+1:(J-W_size+15)*50,:);%SR
% orig_data_M=Mone(15*50+1:(J-W_size+15)*50,:);%AQI
orig_data_M=MAWI(3*50+1:(J-W_size+3)*50,:);%MAWI
% orig_data_M=CAIDA(3*100+1:(J-W_size+3)*100,:);%CAIDA
% orig_data_M=SMD(15*50+1:(J-W_size+15)*50,:);%SMD

% calculate all data's NMAE,Cos for each metric
[p_NMAEs,p_COSes]=getPerformanceNC_orign(orig_data_M, RM_orign);

% calculate NMSE of un-sampled data
unOmega=ones(I,J,K)-omega;
p_unOmega_NMAEs=zeros(1,K);
for k=1:K
    Mk = orig_data_M(:,k);
    RMk = RM_orign(:,k);
    unOmegak(:,:) = unOmega(:,W_size+1:J,k);
    fenmuk = sum(abs(Mk.*unOmegak(:)));
    if fenmuk>0
        p_unOmega_NMAEs(1,k) = sum(abs(RMk.*unOmegak(:)-Mk(:).*unOmegak(:)))/fenmuk;
    else
        p_unOmega_NMAEs(1,k) = 0;
    end
end
             
% calculate total sample ratio
samples_total=sum(sum(sum(omega)));
sampleRatio_total=samples_total/(I*J*K);
        
% calculate sample ratio of real-time monitoring phase
samples_sw=samples_total-W_size*I*K;
sampleRatio_sw=samples_sw/(I*(J-W_size)*K);

%----- finemon -----
%----- input: M,W_size,I,J,K,theta,yita,beta
%----------- M: input tensor,  note thta: the top W_size slices are the result of ST
%----------- W_size: window size
%----------- I J K: M's size
%----------- theta: a uniform constants for constants (c,\mu_0, \delta) in sampling bound
%----------- yita: the threshold of residual(.)
%----------- beta: the parameter in delay-recovery strategy
%----- output: 
%----------- R: recovered data
%----------- omega: sampling positions
%----------- pareto_oemga: sampling positions resulted by refined frequency
% function [R,omega,h_incoms,estimators,rs,ms,com, pareto_oemga]=finemon(M,W_size,I,J,K,theta,yita,beta, isRefine,  epsilon_delta, epsilon_gamma)
function [R,omega,h_incoms,estimators,rs,ers,ms,com, pareto_oemga_large,pareto_oemga_low]=finemon(M,W_size,I,J,K,theta,yita,beta, isRefine,  epsilon_delta, epsilon_gamma,p)
    %initial variables
    W=[]; % window containing the slice index
    U_W_t=[]; % subspace
    U_W_t_indep=[];
    U_union_t=[]; % subspace in delay-recovery strategy
    U_W_index=[];
    U_W_t_indep_index=[];
    rank_W=0; % rank

    R = gpuArray(zeros(I,J,K)); % recovered tensor
    omega=zeros(I,J,K); % sampling positions at input tensor M
    
    Pu=[];
    pesudoinverse_U_W_omega_t=[];
    Pu_omega=[];
    estimators=zeros(1,J);
    com=zeros(1,J);
    count_map=containers.Map;

    % TFA parameters in Section 4.2
    pareto_delta = zeros(K,1); % threshold vector \delta 
    pareto_x_min = zeros(K,1); % minimum vector x_m 
    pareto_alpha = zeros(K,1); % pareto parameter alpha_m
    pareto_L = 0; % extra number of samples in (t,t+inv)
    pareto_oemga=zeros(I,J,K); % the sample positions taken by the refined frequency
    pareto_oemga_large=zeros(I,J,K);
    pareto_oemga_low=zeros(I,J,K);

    rs=zeros(1,J);
    ers=zeros(1,J);
    ms=zeros(1,J);
    times=zeros(1,J);% time cost for each slice (from sampling to recovery)
    times_update_rs=zeros(1,J);
    times_estimatorFull=zeros(1,J);
    times_completion=zeros(1,J);
    times_compute_m=zeros(1,J);
    times_oumega=zeros(1,J);
    times_fullcompletion=zeros(1,J);
    times_fullfor=zeros(1,J);
    times_train2=zeros(1,J);
    incomplete_spaceSlices=zeros(1,J);
    incomplete_spaceSlice=0;
    h_incoms=[];
    beta_count=0;
    Mx=[];
    OptimizationIDX=zeros(1,J);
   
    ischange_U_W=1;
    isSame_UWTindep=0;

    for i=1:J
        %***** initial phase: take full sampling
        if length(W)<W_size
            st_train=tic;
            W=[W i];
            Mi=M(:,i,:);
            R(:,i,:)=Mi;
            com(1,i)=-1;
            omega(:,i,:)=ones(I,K);
            ms(1,i)=I*K;

            if size(U_W_t_indep,2)==0
                estimator=1;
            else
                st_train2=tic;
                if ~isSame_UWTindep
                    Pu_train2=tprod(U_W_t_indep,tpinv(U_W_t_indep));
                end
                estimator_err = Mi-tprod(Pu_train2,Mi);
                estimator=(norm(estimator_err(:,:),'fro')^2)/(norm(Mi(:,:),'fro')^2);
                times_train2(1,i)=toc(st_train2);
            end
            normalize=norm(Mi(:));
            if normalize==0
                tempU=Mi;
            else
                tempU=Mi/normalize;
            end
            if size(U_W_t,2)==0
                n=1;
            else
                n=size(U_W_t,2)+1;
            end
             U_W_t(:,n,:)=tempU;
            if size(U_W_t,2)==0
               U_W_t = gpuArray(U_W_t);
            end
            U_W_index=[U_W_index i];
            estimators(1,i)=estimator;
            if estimator>yita
                U_W_t_indep(:,size(U_W_t_indep,2)+1,:)=tempU;
                U_W_t_indep_index=[U_W_t_indep_index i];
                isSame_UWTindep=0;
                rank_W=rank_W+1;
            else
                isSame_UWTindep=1;% avoid the redundant calculation for Pu_train2
            end
            time=toc(st_train);
            times(1,i)=time;
            rs(1,i)=rank_W;
            ers(1,i)=rank_W;
        else
        
        %***** real-time monitoring phase: slide-window moving, TFA+ESTC
            st_slideW=tic;
            %------ subroutine_sw_removeFromW -----
            [W,U_W_t,U_W_index,ischange_U_W,rank_W] = subroutine_sw_removeFromW(W,U_W_t,U_W_index,count_map,ischange_U_W,rank_W,yita,i);
            times_r=toc(st_slideW);
            times_update_rs(1,i)=times_r;
            ers(1,i)=rank_W;

            if incomplete_spaceSlice>0
                st_full1=tic;
                Mi=M(:,i,:);
                R(:,i,:)=Mi;
                omega(:,i,:)=ones(I,K);
                ms(1,i)=I*K;
                com(:,i)=-1;
                beta_count = beta_count+1;
                Mx(size(Mx,1)+1:size(Mx,1)+I,:,:)=Mi;
                if beta_count==beta
                    beta_count=0;
                    %----- subroutine_sw_STupdateU -----
                    [U_W_t,U_union_t,ischange_U_W,estimators,count] = subroutine_sw_STupdateU(Mx,U_W_t,U_union_t,i,I,ischange_U_W,estimators,yita);
                    
                    time=toc(st_full1);
                    times(1,i)=time;
                    if count>0
                        rank_W=rank_W+count;
                        U_W_index=[U_W_index i];
                        keys=count_map.keys;
                        keys=[keys i];
                        values=count_map.values;
                        values=[values count];
                        count_map=containers.Map(keys,values);
                    end
                    Mx=[];

                    %----- subroutine_delayRecovery -----
                    [R,com,h_incoms,incomplete_spaceSlice,times,times_fullcompletion] = subroutine_delayRecovery(M,R,U_union_t,omega,h_incoms,com,i,incomplete_spaceSlice,times,times_fullcompletion,yita);
           
                    %----- subroutine_delayRecovery_all -----
                    [R,com,h_incoms,times,times_fullfor] = subroutine_delayRecovery_all(M,R,U_union_t,omega,h_incoms,com,i,incomplete_spaceSlice,times,times_fullfor,yita);
                    incomplete_spaceSlice=0;

                    if isempty(h_incoms)
                        U_union_t=[];
                    end
                else
                    time=toc(st_full1);
                    times(1,i)=time;
                end

            else % if incomplete_spaceSlice=0
                st_sample=tic;
                % sampling the ith slice by the frequency
%                 if rank_W==1
%                     m=min(ceil(theta*rank_W*rank_W),I*K);
%                 else
%                     m=min(ceil(theta*rank_W*rank_W*log(rank_W)),I*K);
%                 end
                times_compute_m(1,i)=toc(st_sample);
                m=floor(p*I); %定频采样
                ms(1,i)=m*K;
                omega_column=Get_Array_equalInterval(I,m);
                omega_slice_index=find(omega_column==1);
                % refine the frequency when a large tube is detected
                if isRefine==1 && m~=I
%                     [omega_column, pareto_oemga, L] = refineSampling(omega_slice_index,omega_column, omega, pareto_oemga, M, i, I,K, W,m, epsilon_delta, epsilon_gamma);
                    [omega_column, pareto_oemga_large,pareto_oemga_low, ~] = refineSampling(omega_slice_index,omega_column, omega, pareto_oemga_large,pareto_oemga_low, M, i, I,K, W,m, epsilon_delta, epsilon_gamma);
                
                end
                omega_slice=zeros(I,K);
                for j=1:K
                    omega_slice(:,j)=omega_column;
                end
                % record the sample positions in this slice
                omega(:,i,:)=omega_slice;
                omega_slice_index=find(omega_column==1);
                Mi_omega_1(:,:)=M(:,i,:);
                Mi_omega_1=Mi_omega_1.*omega_slice;
                Mi_omega=[];
                Mi_omega(:,1,:)=Mi_omega_1;
                times_oumega(1,i)=toc(st_sample);
                
                R(:,i,:)=Mi_omega;
                st_estimatorFull=tic;
                % calculate residual(.) = estimator
                isSame_sample=0;
                %----- subroutine_residual -----
                [estimator,pesudoinverse_U_W_omega_t,Pu_omega,isSame_sample,OptimizationIDX] = subroutine_residual(Mi_omega,omega,omega_slice,omega_slice_index,i,U_W_t,pesudoinverse_U_W_omega_t,Pu_omega,isSame_sample,ischange_U_W,W_size,OptimizationIDX);

                times_estimatorFull(1,i)=toc(st_estimatorFull);
                estimators(1,i)=estimator;
                % ESTC: Case Yes if estimator<=yita; Case No if estimator>yita
                if estimator>yita
                    st_estimatorFull2=tic;
                    % Wo TFA：考虑要不要去除全采样，去不去找个说法都行
                    if i==W_size+1 || ~isSame_sample
                        % ridge regression, \lambda = 0.01
                        U_W_omega=U_W_t(omega_slice_index,:,:);
                        ridge = tprod(tran(U_W_omega),U_W_omega);
                        ridge_eye = 0.0001*teye(size(ridge,1), size(ridge,3));
                        pesudoinverse_U_W_omega_t=tprod(tinv(ridge-ridge_eye),tran(U_W_omega)); 

                        Pu=tprod(U_W_t,pesudoinverse_U_W_omega_t);
                        ischange_U_W=1;
                    else
                        OptimizationIDX(1,i)=OptimizationIDX(1,i)+1;
                    end
                    %----- subroutine_ESTC -----
                    [R,com]=subroutine_ESTC(R,Pu,Mi_omega,i,omega_slice_index,com);
                    times_completion(1,i)=toc(st_estimatorFull2);
                    rank_W=rank_W-1;
%                     if m~=I
%                         incomplete_spaceSlice=i;
%                         h_incoms=[h_incoms i];
%                         incomplete_omega_column=omega_column;
%                         beta_count=0;
%     
%                         if size(U_union_t,2)==0
%                             U_union_t = U_W_t;
%                         end
%                     else
%                         normalize=norm(Mi_omega_1(:));
%                         if normalize==0
%                             tempU=Mi_omega_1;
%                         else
%                             tempU=Mi_omega_1./normalize;
%                         end
%                         U_W_t(:,size(U_W_t,2)+1,:)=tempU;
%                         U_W_index=[U_W_index i];
%                         if size(U_union_t,2)>0
%                             U_union_t(:,size(U_union_t,2)+1,:)=tempU;
%                         end
%                     end

                    times_estimatorFull(1,i)=times_estimatorFull(1,i)+toc(st_estimatorFull2);
                else
                    st_completion=tic;

                    if i==W_size+1 || ~isSame_sample
                        % ridge regression, \lambda = 0.01
                        U_W_omega=U_W_t(omega_slice_index,:,:);
                        ridge = tprod(tran(U_W_omega),U_W_omega);
                        ridge_eye = 0.0001*teye(size(ridge,1), size(ridge,3));
                        pesudoinverse_U_W_omega_t=tprod(tinv(ridge-ridge_eye),tran(U_W_omega)); 

                        Pu=tprod(U_W_t,pesudoinverse_U_W_omega_t);
                        ischange_U_W=1;
                    else
                        OptimizationIDX(1,i)=OptimizationIDX(1,i)+1;
                    end
                    %----- subroutine_ESTC -----
                    [R,com]=subroutine_ESTC(R,Pu,Mi_omega,i,omega_slice_index,com);
                    times_completion(1,i)=toc(st_completion);
                    rank_W=rank_W-1;
                end

                time=toc(st_sample);
                times(1,i)=time;
            end %if incomplete_spaceSlice>0
            rs(1,i)=rank_W;
        end
        incomplete_spaceSlices(1,i)=incomplete_spaceSlice;

    end%end for
end