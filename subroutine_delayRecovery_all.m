%----- subroutine_delayRecovery_all: recover the previous all incomplete slices recorded in h_incoms
%----- main inputs: U_union_t, h_incoms
%----- main outputs: R recovered data
function [R,com,h_incoms,times,times_fullfor]=subroutine_delayRecovery_all(M,R,U_union_t,omega,h_incoms,com,i,incomplete_spaceSlice,times,times_fullfor,yita)
    h_incoms2=setdiff(h_incoms,incomplete_spaceSlice);
    for j=1:length(h_incoms2)
        jj=h_incoms2(j);
        st_fullforj=tic;

        isSameOmega=0;
        if j>1
            for jo=1:size(omega(:,jj,:),1)
                if omega(jo,jj,1)==omega(jo,h_incoms2(j-1),1)
                    isSameOmega=1;
                else
                    isSameOmega=0;
                    break
                end
            end
        end
        omega_slice_index=find(omega(:,jj,1)==1);
        isSame_forall=j>1 && isSameOmega;
        %isSame_forall=0;
        if ~isSame_forall
            U_all_omega_t=U_union_t(omega_slice_index,:,:);
            pesudoinverse_U_all_omega_t=tpinv(U_all_omega_t);
            Pu_incom_omega=tprod(U_all_omega_t,pesudoinverse_U_all_omega_t);
        end

        Mj_omega_t=[];
        Mj_omega_t(:,1,:)=M(omega_slice_index,jj,:);
        
%         estimator=(norm(tensor(Mj_omega_t-tprod(Pu_incom_omega,Mj_omega_t)))^2)/(norm(tensor(Mj_omega_t))^2);
        estimator_err = Mj_omega_t-tprod(Pu_incom_omega,Mj_omega_t);
        estimator=(norm(estimator_err(:,:),'fro')^2)/(norm(Mj_omega_t(:,:),'fro')^2);
        if estimator<=yita
            if ~isSame_forall
                Pu_Uall_t=tprod(U_union_t,pesudoinverse_U_all_omega_t);
            end
            R(:,jj,:)=tprod(Pu_Uall_t,Mj_omega_t); 
            R(omega_slice_index,jj,:)=Mj_omega_t; 
            com(:,jj)=1;
            h_incoms=setdiff(h_incoms,jj);
        else
            if ~isSame_forall
                Pu_Uall_t=tprod(U_union_t,pesudoinverse_U_all_omega_t);
            end
            R(:,jj,:)=tprod(Pu_Uall_t,Mj_omega_t); 
            R(omega_slice_index,jj,:)=Mj_omega_t; 
            com(:,jj)=1;
        end
        time=toc(st_fullforj);
        times(1,jj)=times(1,jj)+time;
        times_fullfor(1,i)=times_fullfor(1,i)+time;
    end
end