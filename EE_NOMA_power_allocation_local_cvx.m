function [P_i, EXITFLAG, l] = EE_NOMA_power_allocation_local_cvx(Pt, BW, w, R, R_min_JT_user, gamma, rho, P_fix, kappa, PCM)
L_max=15; % Max number of iterations
isJT = false;

[N_users, N_BSs, ~] = size(gamma);


gamma_ib = NaN(N_users,N_BSs);
for bs = 1:N_BSs
    for i = 1:N_users
        gamma_ib(i,bs) = gamma(i,bs,bs);
    end
end

% Calculates the BS that will transmit to the edge user
[~, bs_edge_user] = max(gamma_ib(N_users,:));

% % ========== Feasiblity test =============
% [P_ib, status] = min_power_global_cvx(Pt, BW, w, R, R_min_JT_user, gamma, isJT); % TODO: Do not take the max power into account?
% 
% % Initial power allocation
% P_ib = Pvec2mat(gamma, isJT, P_ib);
% P_ib_old = P_ib;

P_ib = zeros(N_users,N_BSs);
P_ib_old = zeros(N_users,N_BSs);

% if( status2exitflag(status) < 0) %TODO: Does it still make sense?
%     l = 0;
% else
    for l = 0:L_max
        if(l~=0)
            p_aux = P_ib;
            p_aux(isnan(P_ib)) = 0;
            p_old_aux = P_ib_old;
            p_old_aux(isnan(P_ib_old)) = 0;
            %sum(abs(p_aux(:) - p_old_aux(:)).^2)
            if(sum(abs(p_aux(:) - p_old_aux(:)).^2) < 1e-10 || EXITFLAG<=0)
                break;
            end
            P_ib_old = P_ib;
        end
        
        for bs = 1:N_BSs
            %SINR for non-CoMP case
            % TODO: add a parameter to say which BS transmits to the edge
            % user?
            %[ICI, INUI] = interference(gamma, isJT, P_ib(:));
            [ICI, ~] = interference_v2(gamma, isJT, P_ib);

            R_min_i = R(:,bs);
            if (bs == bs_edge_user)
                R_min_i = [R_min_i; R_min_JT_user];
            end
            users_range = 1:(N_users - (bs~=bs_edge_user));
            [P_i, EXITFLAG, ~] = EE_NOMA_power_allocation_local_single_cell_cvx(Pt, BW, w, R_min_i, gamma_ib(users_range,bs), rho, P_fix, kappa, ICI(users_range,bs));
            if(EXITFLAG<0)
                break;
            end
            %[EE, INUI_i, throughput_single_cell, p_tot_single_cell, user_i_data] = single_cell_EE(gamma_ib(users_range,bs), w, BW, P_i, ICI(users_range,bs), rho, kappa, P_fix);
            % Adds the P_i to the matrix of all BSs powers
            if(length(P_i) ~= N_users)
                P_ib(:,bs) = [P_i; NaN];
            else
                P_ib(:,bs) = P_i;
            end
        end
    end
% end
P_i = P_ib(:);
P_i(isnan(P_i)) = 0;


end