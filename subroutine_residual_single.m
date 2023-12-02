%----- subroutine_residual: calculate the residual(Mi_omega, U_W_omega), determine that which case is performed, Case Yes or Case No
%----- main inputs: sampled slice Mi_omega, subspace U_W_omega
%----- main outputs: the result of residual(.) = estimator
function [estimator_singles,pesudoinverse_U_W_omega_t,Pu_omega,isSame_sample,OptimizationIDX]=subroutine_residual_single(Mi_omega,omega,omega_slice,omega_slice_index,i,U_W,pesudoinverse_U_W_omega_t,Pu_omega,isSame_sample,ischange_U_W,W_size,OptimizationIDX)
    estimator_singles = zeros(1,size(omega,3));
    if size(U_W,2)==0
        estimator_singles(1:size(omega,3)) = 1;
    else
        isSameOmega=0;
        for j=1:size(omega_slice,1)
            if omega_slice(j,1)==omega(j,i-1,1)
                isSameOmega=1;
            else
                isSameOmega=0;
                break;
            end
        end
        isSame_sample = (~ischange_U_W) && isSameOmega;
        %isSame_sample = 0;
        if i==W_size+1 || ~isSame_sample
            if ~isempty(omega_slice_index)
                U_W_omega = U_W(omega_slice_index,:,:);
                Pu_omega = zeros(size(U_W_omega,1),size(U_W_omega,1),size(U_W_omega,3));
                for j = 1:size(U_W,3)
                    Pu_omega(:,:,j) = U_W_omega(:,:,j)* pinv(U_W_omega(:,:,j)); 
                end
            end
        else
            OptimizationIDX(1,i)=OptimizationIDX(1,i)+1;
        end
        if ~isempty(omega_slice_index)
            for j = 1:size(Pu_omega,3)
                Mi_omega_1 = Mi_omega(omega_slice_index,1,j);
                proRes = Pu_omega(:,:,j)*Mi_omega_1;
                estimator_singles(j) = (norm(Mi_omega_1-proRes)^2)/(norm(Mi_omega_1)^2);
            end 
        end
    end
end