function [p_i_tot, p_i_tot_user_side] = user_power_consumption(p_i, i, bs, gamma, rho, kappa, isJT)
    [N_users, N_BSs, ~] = size(gamma);
    N_inner_users = N_users - 1;
    
    % Calculates the number of users per cluster
    isConvPCM = true; % when true, the PCM model for JT will be the same as for conv.
    %TODO: check if the SIC constraint is considered or not.
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


% %TODO: review this function!!!!!
% function [p_i_tot, p_i_tot_user_side] = user_power_consumption(p_i, i, bs, gamma, rho, P_fix, kappa, PCM, isJT)
%     [N_users, N_BSs, ~] = size(gamma);
%     N_inner_users = N_users - 1;
%     
%     
%     % Calculates the number of users per cluster
%     if(~isJT)
%         J_b = N_inner_users*ones(N_BSs,1);
%         J_b(1) = N_users;
%     else
%         J_b = N_users*ones(N_BSs,1);
%     end
%     
%     if(isJT && i==J_b(bs) && bs~=1 && p_i ~= 0) % only takes into account the additional power for coordination if there is power allocated
%         % TODO: add kappa here?
%         p_i_tot = p_i + P_fix/J_b(bs); % Removes the power expenditure at the user side. TODO: Check it. PCM1 and PCM2 p_i*rho may represent a power expenditure at the BS side for the user.
%     else
%         if(PCM=="proposed")
%             p_i_tot = p_i + p_i*rho + kappa*(J_b(bs)-i) + kappa + P_fix/J_b(bs); % TODO: add p_i*rho here and treat kappa*(J_b(bs)-i) + kappa as the user power expenditure and p_i*rho as the BS power expenditure for the user?
%         else
%             p_i_tot = p_i + p_i*rho + P_fix/J_b(bs);
%         end
%     end
%     
%     % User-side power
%     if(isJT && i==J_b(bs) && bs~=1)
%         p_i_tot_user_side = 0; % Removes the power expenditure at the user side. TODO: Check it. PCM1 and PCM2 p_i*rho may represent a power expenditure at the BS side for the user.
%     else
%         if(PCM=="proposed")
%             p_i_tot_user_side = p_i*rho + kappa*(J_b(bs)-i) + kappa; % TODO: add p_i*rho here and treat kappa*(J_b(bs)-i) + kappa as the user power expenditure and p_i*rho as the BS power expenditure for the user?
%         else
%             p_i_tot_user_side = p_i*rho;
%         end
%     end
% end