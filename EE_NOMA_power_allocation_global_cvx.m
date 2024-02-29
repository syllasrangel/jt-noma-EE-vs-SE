function [P_i, EXITFLAG, l] = EE_NOMA_power_allocation_global_cvx(Pt, BW, w, R, R_min_JT_user, gamma, rho, P_fix, kappa, PCM, isJT, SIC_constraint)
L_max=30; % Max number of iterations

[N_users, N_BSs, ~] = size(gamma);


gamma_ib = NaN(N_users,N_BSs);
for bs = 1:N_BSs
    for i = 1:N_users
        gamma_ib(i,bs) = gamma(i,bs,bs);
    end
end


% ========== Feasiblity test =============
[P_i, status] = min_power_global_cvx(Pt, BW, w, R, R_min_JT_user, gamma, isJT, SIC_constraint);

% Initial power allocation
P_i_old = P_i;


if( status2exitflag(status) < 0)
    l = 0;
else
    for l = 0:L_max
        
        if(l~=0)
            if(abs(system_throughput(w, BW, gamma, isJT, P_i_old)/system_power_consumption(P_i_old, gamma, rho, P_fix, kappa, isJT, true) - system_throughput(w, BW, gamma, isJT, P_i)/system_power_consumption(P_i, gamma, rho, P_fix, kappa, isJT, true)) < 0.0001)
                break;
            end
            P_i_old = P_i;
        end
        
        %SINR for non-CoMP case
        [ICI, INUI] = interference(gamma, isJT, P_i);
        P_ib = Pvec2mat(gamma, isJT, P_i);
        
        SINR_ib = P_ib.*gamma_ib./(ICI + INUI + w);
    
        a_ib = SINR_ib./(1+SINR_ib);
        c_ib = log2(1+SINR_ib) - SINR_ib.*log2(SINR_ib)./(1+SINR_ib);
        c_ib(isnan(c_ib))=0;
        
        c1_b = NaN(N_BSs,1);
        if(isJT)
            sum_aux = 0;
            for bs = 1:N_BSs
                sum_aux = sum_aux + P_ib(N_users,bs).*gamma_ib(N_users,bs);
            end
            for bs = 1:N_BSs
                c1_b(bs) = P_ib(N_users,bs).*gamma_ib(N_users,bs)./sum_aux;
            end
        end
        
        [P_i, status, ~] = dinkelbach_algorithm_cvx(Pt, BW, w, R, R_min_JT_user, gamma, rho, P_fix, kappa, PCM, a_ib, c_ib, c1_b, isJT, SIC_constraint);
    end
end
EXITFLAG = status2exitflag(status);

end