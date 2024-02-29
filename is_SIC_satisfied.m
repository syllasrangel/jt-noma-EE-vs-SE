function satisfied = is_SIC_satisfied(P_i, BW, w, R, R_min_JT_user, gamma, isJT)
    [N_users, N_BSs, ~] = size(gamma);
    N_inner_users = N_users - 1;
    
    
    % Calculates the number of users per cluster
    P_ib = Pvec2mat(gamma, isJT, P_i);
    J_b = zeros(N_BSs,1);
    for b = 1:N_BSs
       J_b(b) = sum(P_ib(:,b)>10e-9);
    end
%     if(~isJT)
%         J_b = N_inner_users*ones(N_BSs,1);
%         J_b(1) = N_users;
%     else
%         J_b = N_users*ones(N_BSs,1);
%     end
    
    INUI = zeros(N_inner_users+1, N_BSs);
    ICI = zeros(N_inner_users+1, N_BSs);
    INUI_jk = zeros(N_inner_users+1, N_BSs, N_inner_users+1);
    
    [ICI, INUI, INUI_jk] = interference_CVX(gamma, isJT, P_i, ICI, INUI, INUI_jk);
    
    satisfied = true;
    for bs = 1:N_BSs
        for j = 1:(J_b(bs)-1)
            for k = j+1:J_b(bs)
    %     for j = 1:(J_b(bs)-1-isJT) % Not testing conditions to decode edge user`s signal in any user (it allows edge user to be served by only one BS)
    %         for k = j+1:J_b(bs)-isJT % Not testing conditions to decode edge user`s signal in any user (it allows edge user to be served by only one BS)
                if((isJT && k==J_b(bs)) || (~isJT && k==J_b(bs) && bs==1))
                    gamma_min = 2^(R_min_JT_user/(w*BW)) - 1;
                else
                    gamma_min = 2^(R(k,bs)/(w*BW)) - 1;
                end
                %-(q_test(two_dim_2_one_dim(k, bs, N_users, isJT)) - log(ICI(j,bs) + INUI_jk(j,bs,k) + w)./log(2) + log(gamma(j,bs,bs)/gamma_min)./log(2)) <= 0;
                satisfied = satisfied && P_i(two_dim_2_one_dim(k, bs, N_users, isJT))*gamma(j,bs,bs) - gamma_min*(ICI(j,bs) + INUI_jk(j,bs,k) + w) >= 0;
            end
        end
    end
end