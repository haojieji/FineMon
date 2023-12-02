%refine the frequency to capture extreme tubes
% function [omega_column, pareto_oemga, L] = refineSampling(omega_slice_index,omega_column, omega, pareto_oemga, M, i, I,K,W,m, epsilon_delta, epsilon_gamma)
function [omega_column, pareto_oemga_large, pareto_oemga_low, L] = refineSampling(omega_slice_index,omega_column, omega, pareto_oemga_large, pareto_oemga_low, M, i, I,K,W,m, epsilon_delta, epsilon_gamma)
    for j=1:length(omega_slice_index)
        if omega_slice_index(j)<I
            omega(omega_slice_index(j),i,:)=1;
            pareto_delta = zeros(1,K);
            pareto_gamma = zeros(1,K);
            pareto_alpha = zeros(1,K);
            pareto_x_min = zeros(1,K);

            for k = 1:K
                MWk = M(:,W,k);
                omegaWk = omega(:,W,k);
                mWk = MWk(:);
                oWk = omegaWk(:);
                oWk_IDX = find(oWk == 1);
                %pareto_delta(k)=2.7*mean(mWk(oWk_IDX));%CAIDA
                %pareto_delta(k)=1.8*mean(mWk(oWk_IDX));%MAWI
                pareto_delta(k)=epsilon_delta * mean(mWk(oWk_IDX));
                pareto_gamma(k)=epsilon_gamma * mean(mWk(oWk_IDX));
    
                mWk_oWk = mWk(oWk_IDX);
                
                pareto_x_min(k) = min(mWk_oWk(mWk_oWk>0));
                pareto_alpha(k) = length(oWk_IDX)/(sum( log( mWk(oWk_IDX)/pareto_x_min(k) ) ));
            end
            m_ji = M(omega_slice_index(j),i,:);
            is_large = m_ji(:,:) > pareto_delta;
            is_low = m_ji(:,:) < pareto_gamma;
    
            L = 0;
            if max(is_large)==1 && max(is_low)==1
                % contain both large and low metric, we calculete L by the equation (16) in Section 4.2
                min_pareto_delta = min(pareto_delta);
                pareto_x_min_delta = pareto_x_min/min_pareto_delta;
                pareto_x_min_delta_alpha = power(pareto_x_min_delta,pareto_alpha);
                L_large =  floor( min(I-omega_slice_index(j), (ceil(I/m)-1))*max(pareto_x_min_delta_alpha) );
    
                max_pareto_gamma = max(pareto_gamma);
                pareto_x_min_gamma = pareto_x_min/max_pareto_gamma;
                pareto_x_min_gamma_alpha = 1 - power(pareto_x_min_gamma, pareto_alpha);
                L_low =  floor( min(I-omega_slice_index(j), (ceil(I/m)-1))*max(pareto_x_min_gamma_alpha) );
                
                L = max(L_large, L_low);
            elseif max(is_large)==1
                min_pareto_delta = min(pareto_delta);
                pareto_x_min_delta = pareto_x_min/min_pareto_delta;
                pareto_x_min_delta_alpha = power(pareto_x_min_delta,pareto_alpha);
                L =  floor( min(I-omega_slice_index(j), (ceil(I/m)-1))*max(pareto_x_min_delta_alpha) );
            elseif max(is_low)==1
                max_pareto_gamma = max(pareto_gamma);
                pareto_x_min_gamma = pareto_x_min/max_pareto_gamma;
                pareto_x_min_gamma_alpha = 1 - power(pareto_x_min_gamma, pareto_alpha);
                L =  floor( min(I-omega_slice_index(j), (ceil(I/m)-1))*max(pareto_x_min_gamma_alpha) );
            end
    
            % sampling extra L tubes in interval (omega_slice_index(j), omega_slice_index(j+1) or I)
            if L~=0
                Ll = min(I-omega_slice_index(j), (ceil(I/m)-1));
                if Ll<L
                    L=Ll;
                end
                omega_column_L = Get_Array_equalInterval(Ll,L);
                omega_column( omega_slice_index(j)+1 : omega_slice_index(j)+Ll ) = omega_column_L;
                omega( omega_slice_index(j) + find(omega_column_L==1), i, :) = 1;
%                 pareto_oemga(omega_slice_index(j)+find(omega_column_L==1),i,:) = 1;
                if max(is_large)==1 && max(is_low)==1
                    if L_large>L_low
                        pareto_oemga_large(omega_slice_index(j)+find(omega_column_L==1),i,:) = 1;
                    else
                        pareto_oemga_low(omega_slice_index(j)+find(omega_column_L==1),i,:) = 1;
                    end
                elseif max(is_large)==1
                    pareto_oemga_large(omega_slice_index(j)+find(omega_column_L==1),i,:) = 1;
                elseif max(is_low)==1
                    pareto_oemga_low(omega_slice_index(j)+find(omega_column_L==1),i,:) = 1;
                end
            end
        end
    end