%----- paper_subroutine_delayRecovery: recover the overdated incomplete_spaceSlice-th slice using the union subspace U_union_t
%----- main inputs: U_union_t, incomplete_spaceSlice
%----- main outputs: the recovered data R
function [R,com,h_incoms,incomplete_spaceSlice,times,times_fullcompletion]=subroutine_delayRecovery(M,R,U_union_t,omega,h_incoms,com,i,incomplete_spaceSlice,times,times_fullcompletion,yita)
                
    incomplete_spaceSlice_index=find(omega(:,incomplete_spaceSlice,1)==1);
    U_incom_omega_t = U_union_t(incomplete_spaceSlice_index,:,:);
    st_fullcompletion=tic;
    pesudoinverse_U_incom_omega_t=tpinv(U_incom_omega_t);
    Pu_incom_omega=tprod(U_incom_omega_t,pesudoinverse_U_incom_omega_t);

    Mi_omega_t=M(incomplete_spaceSlice_index,incomplete_spaceSlice,:);
%     estimator=(norm(tensor(Mi_omega_t-tprod(Pu_incom_omega,Mi_omega_t)))^2)/(norm(tensor(Mi_omega_t))^2);
    estimator_err = Mi_omega_t-tprod(Pu_incom_omega,Mi_omega_t);
    estimator=(norm(estimator_err(:,:),'fro')^2)/(norm(Mi_omega_t(:,:),'fro')^2);
    if estimator<=yita
        alphai=tprod(pesudoinverse_U_incom_omega_t,Mi_omega_t);
        R(:,incomplete_spaceSlice,:)=tprod(U_union_t,alphai);     
        R(incomplete_spaceSlice_index,incomplete_spaceSlice,:)=Mi_omega_t;
        com(:,incomplete_spaceSlice)=1;

        time=toc(st_fullcompletion);
        times(1,incomplete_spaceSlice)=times(1,incomplete_spaceSlice)+time;
        times_fullcompletion(1,i)=time;
        h_incoms=setdiff(h_incoms,incomplete_spaceSlice);
        incomplete_spaceSlice=0;
    else
        alphai=tprod(pesudoinverse_U_incom_omega_t,Mi_omega_t);
        R(:,incomplete_spaceSlice,:)=tprod(U_union_t,alphai);     
        R(incomplete_spaceSlice_index,incomplete_spaceSlice,:)=Mi_omega_t;
        com(:,incomplete_spaceSlice)=1;

        time=toc(st_fullcompletion);
        times(1,incomplete_spaceSlice)=times(1,incomplete_spaceSlice)+time;
        times_fullcompletion(1,i)=time;
    end
end
