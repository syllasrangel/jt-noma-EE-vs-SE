function [p_i_tot, p_i_tot_user_side] = user_power_consumption(p_i, i, bs, gamma, rho, kappa, isJT)
    [N_users, N_BSs, ~] = size(gamma);
    N_inner_users = N_users - 1;
    
    % Calculates the number of users per cluster
    isConvPCM = true; % when true, the PCM model for JT will be the same as for conv.
    if(isConvPCM || ~isJT)
        J_b = N_inner_users*ones(N_BSs,1);
        J_b(1) = N_users;
    else
        J_b = N_users*ones(N_BSs,1);
    end
    
    % For convenience, in JT, the decoding power is added only at the
    % (user = J_b(bs), BS = 1) pair. This is done to avoid accounting more
    % than once for the decoding power, since (user = J_b(bs), BS > 1) 
    % represents the same user in JT and this users will always decode only
    % once.
    if(isJT && i==J_b(bs)+isConvPCM && bs~=1) 
        p_i_tot = p_i + p_i*rho;
    else
        p_i_tot = p_i + p_i*rho + kappa*(J_b(bs)-i+1); 
    end
    
    
    % User-side power
    if(isJT && i==J_b(bs)+isConvPCM && bs~=1)
        p_i_tot_user_side = p_i*rho;
%     elseif(isJT && i==J_b(bs) && bs==1) 
%         p_i_tot_user_side = p_i*rho + kappa;
    else
        p_i_tot_user_side = p_i*rho + kappa*(J_b(bs)-i+1); 
    end
end