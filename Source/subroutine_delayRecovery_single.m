%----- paper_subroutine_delayRecovery: recover the overdated incomplete_spaceSlice-th slice using the union subspace U_union_t
%----- main inputs: U_union_t, incomplete_spaceSlice
%----- main outputs: the recovered data R
function [R,com,h_incoms,incomplete_spaceSlice,times,times_fullcompletion]=subroutine_delayRecovery_single(M,R,U_union_t,omega,h_incoms,com,i,incomplete_spaceSlice,times,times_fullcompletion,yita)
                
    incomplete_spaceSlice_index=find(omega(:,incomplete_spaceSlice,1)==1);
    U_incom_omega_t = U_union_t(incomplete_spaceSlice_index,:,:);
    st_fullcompletion=tic;

    Pu_incom_omega = zeros(size(U_incom_omega_t,1), size(U_incom_omega_t,1), size(U_incom_omega_t,3));
    for j=1:size(U_union_t,3)
        Pu_incom_omega(:,:,j) = U_incom_omega_t(:,:,j)*pinv(U_incom_omega_t(:,:,j));
    end

    Mi_omega_t=M(incomplete_spaceSlice_index,incomplete_spaceSlice,:);
    for j=1:size(U_union_t,3)
        estimator=(norm(Mi_omega_t(:,:,j)-Pu_incom_omega(:,:,j)*Mi_omega_t(:,:,j))^2)/(norm(Mi_omega_t(:,:,j))^2);

        if estimator<=yita
            R(:,incomplete_spaceSlice,j) = U_union_t(:,:,j)*pinv(U_incom_omega_t(:,:,j))*Mi_omega_t(:,1,j);
            R(incomplete_spaceSlice_index,incomplete_spaceSlice,j) = Mi_omega_t(:,:,j);
            com(:,incomplete_spaceSlice)=1;

            h_incoms=setdiff(h_incoms,incomplete_spaceSlice);
            
        else
            R(:,incomplete_spaceSlice,j) = U_union_t(:,:,j)*pinv(U_incom_omega_t(:,:,j))*Mi_omega_t(:,1,j);
            R(incomplete_spaceSlice_index,incomplete_spaceSlice,j) = Mi_omega_t(:,1,j);
            com(:,incomplete_spaceSlice)=1; 
        end
    end
    time=toc(st_fullcompletion);
    times(1,incomplete_spaceSlice)=times(1,incomplete_spaceSlice)+time;
    times_fullcompletion(1,i)=time;
    if ~ismember(h_incoms, incomplete_spaceSlice)
        incomplete_spaceSlice=0;
    end
end
