function [P_i, status] = min_power_global_cvx(Pt, BW, w, R, R_min_JT_user, gamma, isJT)
[N_users, N_BSs, ~] = size(gamma);
N_inner_users = N_users - 1;

obj_func = @(x)sum(x);

length_Pi = N_BSs*(N_inner_users + isJT) + ~isJT;

% Calculates the number of users per cluster
if(~isJT)
    J_b = N_inner_users*ones(N_BSs,1);
    J_b(1) = N_users;
else
    J_b = N_users*ones(N_BSs,1);
end


cvx_begin quiet
    variable P_i(length_Pi);
    expression INUI(N_inner_users+1, N_BSs)
    expression ICI(N_inner_users+1, N_BSs)
    minimize( obj_func(P_i) );
    subject to
        % ICI and inter-NOMA-user interference (INUI) calculation
        [ICI, INUI] = interference_CVX(gamma, isJT, P_i, ICI, INUI);

        % === Rate requirement inner users ===
        for bs = 1:N_BSs
            for j = 1:N_inner_users
                gamma_min = 2^(R(j,bs)/(w*BW)) - 1;
                P_i(two_dim_2_one_dim(j, bs, N_users, isJT))*gamma(j,bs,bs) - gamma_min*(ICI(j,bs) + INUI(j,bs) + w) >= 0;
            end
        end

        % === Rate requirement cell-edge user ===
        gamma_min = 2^(R_min_JT_user/(w*BW)) - 1;
        if(~isJT)
            P_i(two_dim_2_one_dim(J_b(1), 1, N_users, isJT))*gamma(J_b(1),1,1) - gamma_min*(ICI(J_b(1),1) + INUI(J_b(1),1) + w) >= 0;
        else
            sum_useful_term = 0;
            for bs_aux = 1:N_BSs % Relaxed constraint for JT case.
                sum_useful_term = sum_useful_term + P_i(two_dim_2_one_dim(J_b(bs_aux), bs_aux, N_users, isJT))*gamma(J_b(bs_aux),bs_aux,bs_aux);
            end 
            sum_useful_term - gamma_min*(ICI(J_b(1),1) + INUI(J_b(1),1) + w) >= 0;
        end
         % === BS max power constraint ===
        for bs = 1:N_BSs
            tot_power = 0;
            for j = 1:J_b(bs) 
                tot_power = tot_power + P_i(two_dim_2_one_dim(j, bs, N_users, isJT));
            end
            tot_power - Pt <= 0;
        end
        % === Non-negative power ===
        P_i >= 0;
cvx_end

% cvx_status
% cvx_optval

% Sometimes a small negative value (arround -1e-14) is possible due to the
% solver's tolerance (maeke it zero!)
P_i(P_i < 0) = 0;

status = cvx_status;

end