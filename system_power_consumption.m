function Pt = system_power_consumption(P_i, gamma, rho, P_fix, kappa, isJT, is_cvx)
    [N_users, N_BSs, ~] = size(gamma);
    N_inner_users = N_users - 1;
    
    force_PCM_conv = true; % forces the PCM conv in the cvx allocation
%     % Calculates the number of users per cluster
%     if(~isJT)
%         J_b = N_inner_users*ones(N_BSs,1);
%         J_b(1) = N_users;
%     else
%         J_b = N_users*ones(N_BSs,1);
%     end
%     
%     two2one = @(j,bs)two_dim_2_one_dim(j, bs, N_users, isJT);
       
    %P_ib = Pvec2mat(gamma, isJT, P_i);
    P_ib = Pvec2mat(gamma, length(P_i) == N_BSs*N_users, P_i);
    
%     Pt = 0;
    Pt = 0;
    for bs = 1:N_BSs
%         Pt = Pt + P_fix;
        Pt = Pt + P_fix + sum(user_power_consumption_2(P_ib(:,bs), gamma, rho, kappa, is_cvx));
%         for i = 1:J_b(bs)
%             Pt = Pt + user_power_consumption(P_i(two2one(i,bs)), i, bs, gamma, rho, kappa, isJT);
%         end
    end
    % Workaround: for PCM proposed in CVX. 
    % If isJT is true and force_PCM_conv false, it will calculate the SIC
    % overhead for edge user in all BS even if the allocated power is zero
    % in some BS.
    if(is_cvx && ~force_PCM_conv && isJT)
        Pt = Pt + N_BSs*N_users*kappa;
    elseif(is_cvx && (~isJT || force_PCM_conv)) % if is CVX and conv. then only one BS transmit for the edge user.
        Pt = Pt + N_users*kappa;
    end
end