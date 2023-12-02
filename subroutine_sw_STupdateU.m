%----- subroutine_sw_STupdateU: take self-embedding transform for full sampling slices (Section 4.3.1), update subspace U_W_t
function [U_W_t,U_union_t,ischange_U_W,estimators,count] = subroutine_sw_STupdateU(Mx,U_W_t,U_union_t,i,I,ischange_U_W,estimators,yita)
    count=0;
    isSameUw_full=0;
    isSameUall_full=0;
    for j=1:size(Mx,1)-I+1
        MXj=Mx(j:j+I-1,:,:);
        % update U_W
        if size(U_W_t,2)==0
            estimator=1;
        else
            if ~isSameUw_full
                Pu_full=tprod(U_W_t,tpinv(U_W_t));
            end
            estimator=(norm(tensor(MXj-tprod(Pu_full,MXj)))^2)/(norm(tensor(MXj))^2);
        end
        if estimator>yita
            normalize=norm(tensor(MXj));
            if normalize==0
                tempU=MXj;
            else
                tempU=MXj/normalize;
            end
            if size(U_W_t,2)==0
                n=1;
            else
                n=size(U_W_t,2)+1;
            end
            U_W_t(:,n,:)=tempU;
            ischange_U_W=1;
            isSameUw_full=0;
            estimators(1,i)=estimator;
            count=count+1;
        else
            isSameUw_full=1;
        end

        if ~isSameUall_full
            Pu_all_full=tprod(U_union_t,tpinv(U_union_t));
        end
        estimator=(norm(tensor(MXj-tprod(Pu_all_full,MXj)))^2)/(norm(tensor(MXj))^2);
        if estimator>yita
            normalize=norm(tensor(MXj));
            if normalize==0
                tempU=MXj;
            else
                tempU=MXj/normalize;
            end
            if size(U_union_t,1)==0
                n=1;
            else
                n=size(U_union_t,2)+1;
            end
            U_union_t(:,n,:)=tempU;
            isSameUall_full=0;
        else
            isSameUall_full=1;
        end
    end
end