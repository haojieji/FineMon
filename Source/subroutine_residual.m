%----- subroutine_residual: calculate the residual(Mi_omega, U_W_omega), determine that which case is performed, Case Yes or Case No
%----- main inputs: sampled slice Mi_omega, subspace U_W_omega
%----- main outputs: the result of residual(.) = estimator
function [estimator,pesudoinverse_U_W_omega_t,Pu_omega,isSame_sample,OptimizationIDX]=subroutine_residual(Mi_omega,omega,omega_slice,omega_slice_index,i,U_W,pesudoinverse_U_W_omega_t,Pu_omega,isSame_sample,ischange_U_W,W_size,OptimizationIDX)
    if size(U_W,2)==0
        estimator=1;
    else
        isSameOmega=0;
        for j=1:size(omega_slice,1)
            if omega_slice(j,1)==omega(j,i-1,1)
                isSameOmega=1;
            else
                isSameOmega=0;
                break
            end
        end
        isSame_sample=(~ischange_U_W) && isSameOmega;
        %isSame_sample=0;
        if i==W_size+1 || ~isSame_sample
            if ~isempty(omega_slice_index)
                U_W_omega=U_W(omega_slice_index,:,:);
                pesudoinverse_U_W_omega_t=tpinv(U_W_omega);
                Pu_omega=tprod(U_W_omega,pesudoinverse_U_W_omega_t);
            end
        else
            OptimizationIDX(1,i)=OptimizationIDX(1,i)+1;
        end
        proRes=[];
        if ~isempty(omega_slice_index)
            Mi_omega_1=Mi_omega(omega_slice_index,:,:);
            proRes(:,:)=tprod(Pu_omega,Mi_omega_1);
            estimator=(norm(Mi_omega_1(:,:)-proRes,'fro')^2)/(norm(Mi_omega_1(:,:),'fro')^2);
        else
            estimator=0;
        end
    end
end