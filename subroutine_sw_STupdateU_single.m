%----- subroutine_sw_STupdateU: take self-embedding transform for full sampling slices (Section 4.3.1), update subspace U_W_t
function [U_W_t,U_union_t,ischange_U_W,estimators,count,rank_Ws] = subroutine_sw_STupdateU_single(Mx,U_W_t,U_union_t,i,I,ischange_U_W,estimators,yita,rank_Ws)
    count=0;
    isSameUw_full=0;
    isSameUall_full=0;
    for j=1:size(Mx,1)-I+1
        MXj=Mx(j:j+I-1,:,:);
        % update U_W
        estimator_singles = zeros(1,size(Mx,3));
        if size(U_W_t,2)==0
            estimator_singles = 1;
        else
            if ~isSameUw_full
                Pu_full= zeros(size(U_W_t,1),size(U_W_t,1),size(U_W_t,3));
                for jj = 1:size(U_W_t,3)
                    Pu_full(:,:,jj) = U_W_t(:,:,jj)*pinv(U_W_t(:,:,jj));
                end
            end
            for jj = 1:size(U_W_t,3)
                estimator_singles(jj) = (norm(MXj(:,1,jj)-Pu_full(:,:,jj)*MXj(:,1,jj))^2)/(norm(MXj(:,1,jj))^2);
            end
        end
        updateIndicator = estimator_singles>yita;
        rank_Ws(updateIndicator) = rank_Ws(updateIndicator) +1;
        if ~isempty(find(updateIndicator, 1))
            tempU = zeros(size(U_W_t,1), 1, size(U_W_t,3));
            for jj=1:size(U_W_t,3) 
                normalize = norm(MXj(:,1,jj));
                if normalize==0
                    tempU(:,1,jj) = MXj(:,1,jj);
                else
                    tempU(:,1,jj) = MXj(:,1,jj)./normalize;
                end
            end
            if size(U_W_t,2)==0
                n=1;
            else
                n=size(U_W_t,2)+1;
            end
            U_W_t(:,n,:)=tempU;
            ischange_U_W=1;
            isSameUw_full=0;
            estimators(1,i)=max(estimator_singles);
            count=count+1;
 
        else
            isSameUw_full=1;
        end

        if ~isSameUall_full
            Pu_all_full=zeros(size(U_union_t,1),size(U_union_t,1),size(U_union_t,3));
            for jj = 1:size(U_union_t,3)
                Pu_all_full(:,:,jj) = U_union_t(:,:,jj)*pinv(U_union_t(:,:,jj));
            end
        end
        for jj = 1:size(U_union_t,3)
            estimator_singles(jj) = (norm(MXj(:,1,jj)-Pu_all_full(:,:,jj)*MXj(:,1,jj))^2)/(norm(MXj(:,1,jj))^2);
        end
        if ~isempty(find(estimator_singles>yita, 1))
            tempU = zeros(size(U_W_t,1), 1, size(U_W_t,3));
            for jj=1:size(U_W_t,3) 
                normalize = norm(MXj(:,1,jj));
                if normalize==0
                    tempU(:,1,jj) = MXj(:,1,jj);
                else
                    tempU(:,1,jj) = MXj(:,1,jj)./normalize;
                end
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