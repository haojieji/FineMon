 %------ subroutine_sw_removeFromW: removing the oldest slice from the ative window, updating the rank of common window r, estimating active window's rank r'=r+1
 %------ inputs: W, U_W, U_W_t, U_W_index, count_map, ischange_U_W, rank_W
 %------ outputs: W, U_W, U_W_t, U_W_index, ischange_U_W, rank_W
function [W,U_W_t,U_W_index,ischange_U_W,rank_Ws]=subroutine_sw_removeFromW_single(W,U_W_t,U_W_index,count_map,ischange_U_W,rank_Ws,yita,i)
    oldest_slice_index=W(1);
    W=W(1,2:length(W));
    if ismember(oldest_slice_index,U_W_index)
        if isKey(count_map,oldest_slice_index)
            count=count_map(oldest_slice_index);
            if size(U_W_t,2)>count
                U_W_remove=U_W_t(:,1:count,:);
                U_W_t=U_W_t(:,count+1:size(U_W_t,2),:);
                
                for j=1:size(U_W_t, 3)
                    Pu_slideW2_j = U_W_t(:,:,j)*pinv(U_W_t(:,:,j));
                    for jj = 1:count
                        proRes = Pu_slideW2_j*U_W_remove(:,jj,j);
                        estimator = (norm(U_W_remove(:,jj,j)-proRes)^2)/(norm(U_W_remove(:,jj,j))^2);
                        if estimator>yita
                            ischange_U_W=1;
                            if rank_Ws(j)>1
                                rank_Ws(j)=rank_Ws(j)-1;
                            end
                        end
                    end
                end

                U_W_index=setdiff(U_W_index,oldest_slice_index);
            end
        else
            U_W_remove=[];
            if size(U_W_t,2)>1
                U_W_remove(:,1,:)=U_W_t(:,1,:);
                U_W_t=U_W_t(:,2:size(U_W_t,2),:);
                for j=1:size(U_W_t, 3)
                    Pu_slideW2_j = U_W_t(:,:,j)*pinv(U_W_t(:,:,j));
                    proRes = Pu_slideW2_j*U_W_remove(:,1,j);
                    estimator = (norm(U_W_remove(:,1,j)-proRes)^2)/(norm(U_W_remove(:,1,j))^2);
                    if estimator>yita
                        ischange_U_W=1;
                        if rank_Ws(j)>1
                              rank_Ws(j)=rank_Ws(j)-1;
                        end
                    end
                end

                U_W_index=setdiff(U_W_index,oldest_slice_index);
            end
        end
    end

    rank_Ws = rank_Ws + 1;
    W=[W i];
end