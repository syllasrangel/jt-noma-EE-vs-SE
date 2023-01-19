% q => (N_inner_users+1)*N_BSs x 1 (if JT) and (N_inner_users*N_BSs +1) x 1
% (if conventional).
% q_ib => N_inner_users+1 x N_BSs
% P_ib = 2.^q_ib
function nonlcon_ICI_global_CVX(w, BW, gamma, Ri, R_min_JT_user, Pt, P_tol, isJT, initial, c_b, ICI, INUI, q)
    [N_users, N_BSs, ~] = size(gamma);
    N_inner_users = N_users - 1;
    
    % Calculates the number of users per cluster
    if(~isJT)
        J_b = N_inner_users*ones(N_BSs,1);
        J_b(1) = N_users;
    else
        J_b = N_users*ones(N_BSs,1);
    end
    
    % Reshapes the vector q   
    q_ib = Pvec2mat(gamma, isJT, q);

    % ICI and inter-NOMA-user interference (INUI) calculation
    [ICI, INUI] = interference_CVX(gamma, isJT, 2.^q, ICI, INUI);
    
    idx_counter = 0;
    
    % === SIC constraint ===
%     for bs = 1:N_BSs
%         for k = 1:J_b(bs)-1
%             for ii = k+1:J_b(bs)
%                 sum_INUI_power = 0;
%                 for aa = 1:ii-1
%                     sum_INUI_power = sum_INUI_power + 2.^q_ib(aa,bs).*gamma(k,bs,bs);
%                 end
%                 if(isJT && ii == J_b(bs))
%                     sum_useful_power = 0;
%                     for bs_aux = 1:N_BSs
%                         sum_useful_power = sum_useful_power + 2.^q_ib(ii,bs_aux).*gamma(k,bs_aux,bs_aux);
%                     end 
%                 else
%                     sum_useful_power = 2.^q_ib(ii,bs).*gamma(k,bs,bs);
%                 end
%                 idx_counter = idx_counter+1;
%                 %c(idx_counter) = -(sum_useful_power - sum_INUI_power - ICI(k,bs) - P_tol);
%                 -(sum_useful_power - sum_INUI_power - ICI(k,bs) - P_tol) <= 0;
%             end
%         end
%     end
    
    % === Rate requirement inner users ===
    for bs = 1:N_BSs
        for j = 1:N_inner_users
            gamma_min = 2^(Ri(j,bs)/(w*BW)) - 1;
            idx_counter = idx_counter+1;
            %c(idx_counter) = -(q_ib(j,bs) - log2(ICI(j,bs) + INUI(j,bs) + w) + log2(gamma(j,bs,bs)/gamma_min));
            -(q(two_dim_2_one_dim(j, bs, N_users)) - log(ICI(j,bs) + INUI(j,bs) + w)./log(2) + log(gamma(j,bs,bs)/gamma_min)./log(2)) <= 0;
        end
    end
    
    % === Rate requirement cell-edge user ===
    gamma_min = 2^(R_min_JT_user/(w*BW)) - 1;
    idx_counter = idx_counter+1;
    if(~isJT)
        %c(idx_counter) = -(q_ib(J_b(1),1) - log2(ICI(J_b(1),1) + INUI(J_b(1),1) + w) + log2(gamma(J_b(1),1,1)/gamma_min));
        -(q(two_dim_2_one_dim(J_b(1),1)) - log2(ICI(J_b(1),1) + INUI(J_b(1),1) + w) + log2(gamma(J_b(1),1,1)/gamma_min))<=0;
    else
        sum_useful_term = 0;
        for bs_aux = 1:N_BSs % Relaxed constraint for JT case.
            sum_useful_term = sum_useful_term + q(two_dim_2_one_dim(J_b(bs_aux),bs_aux)).*c_b(bs_aux) + c_b(bs_aux).*log2(gamma(J_b(bs_aux),bs_aux,bs_aux)/c_b(bs_aux));
        end 
        %c(idx_counter) = -(sum_useful_term - log2(ICI(J_b(1),1) + INUI(J_b(1),1) + w) - log2(gamma_min));
        -(sum_useful_term - log2(ICI(J_b(1),1) + INUI(J_b(1),1) + w) - log2(gamma_min))<=0;
    end
    
    % === BS max power constraint ===
    if (~initial)
        for bs = 1:N_BSs
            tot_power = 0;
            for j = 1:N_BSs 
                tot_power = tot_power + 2.^(q(two_dim_2_one_dim(j,bs)));
            end
            idx_counter = idx_counter+1;
            %c(idx_counter) = tot_power - Pt;
            tot_power - Pt <= 0;
        end
    end
      
end