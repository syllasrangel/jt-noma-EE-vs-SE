function [P_i, EXITFLAG, l] = EE_NOMA_power_allocation_global_fmincon(Pt, P_tol, BW, w, R, R_min_JT_user, gamma, rho, P_fix, kappa, PCM, isJT)
L_max=30; % iteration number
%isJT = true;

[N_users, N_BSs, ~] = size(gamma);
N_inner_users = N_users - 1;


A=[];
for bs = 1:N_BSs
   aux = (bs ~= 1 && ~isJT);
   A = [A, [zeros(bs-1,N_users - aux); ones(1, N_users - aux); zeros(N_BSs-bs,N_users - aux)]];
end
b = Pt*ones(N_BSs,1);


[P_i, EXITFLAG1] = min_power_global_fmincon(Pt, BW, w, R, R_min_JT_user, gamma, rho, P_fix, kappa, PCM, isJT);
%[P_i2, status] = min_power_global(Pt, BW, w, R, R_min_JT_user, gamma, rho, P_fix, kappa, PCM, isJT);



gamma_ib = NaN(N_users,N_BSs);
for bs = 1:N_BSs
    for i = 1:N_users
        gamma_ib(i,bs) = gamma(i,bs,bs);
    end
end

% Calculates the number of users per cluster
if(~isJT)
    J_b = N_inner_users*ones(N_BSs,1);
    J_b(1) = N_users;
else
    J_b = N_users*ones(N_BSs,1);
end

%q_test = log2(P_i);
P_i_old = P_i;

if( sum(A*(P_i)<=b) ~= length(b) )
    EXITFLAG = -6;
    l=0;
else
    for l = 0:L_max
        
        if(l~=0)
            %abs(system_throughput(w, BW, gamma, isJT, 2.^q_test)/system_power_consumption(2.^q_test, gamma, rho, P_fix, kappa, PCM, isJT) - system_throughput(w, BW, gamma, isJT, P_i)/system_power_consumption(P_i, gamma, rho, P_fix, kappa, PCM, isJT))
            if(abs(system_throughput(w, BW, gamma, isJT, P_i_old)/system_power_consumption(P_i_old, gamma, rho, P_fix, kappa, PCM, isJT) - system_throughput(w, BW, gamma, isJT, P_i)/system_power_consumption(P_i, gamma, rho, P_fix, kappa, PCM, isJT)) < 0.0001)
                break;
            end
            P_i_old = P_i;
        end
        %P_i = 2.^q_test;
        
        %SINR for non-CoMP case
        [ICI, INUI] = interference(gamma, isJT, P_i);
        P_ib = Pvec2mat(gamma, isJT, P_i);
        
        SINR_ib = P_ib.*gamma_ib./(ICI + INUI + w);
    
        a_ib = SINR_ib./(1+SINR_ib);
        c_ib = log2(1+SINR_ib) - SINR_ib.*log2(SINR_ib)./(1+SINR_ib);
        %c_ib = log2(1+SINR_ib) - log2(SINR_ib.^SINR_ib)./(1+SINR_ib);
        c_ib(isnan(c_ib))=0;
        
        c1_b = NaN(N_BSs,1);
        sum_aux = 0;
        for bs = 1:N_BSs
            sum_aux = sum_aux + P_ib(N_users,bs).*gamma_ib(N_users,bs);
        end
        for bs = 1:N_BSs
            c1_b(bs) = P_ib(N_users,bs).*gamma_ib(N_users,bs)./sum_aux;
        end
        
        [P_i, EXITFLAG, l_dinkelbach] = dinkelbach_algorithm_fmincon(Pt, BW, w, R, R_min_JT_user, gamma, rho, P_fix, kappa, PCM, a_ib, c_ib, c1_b, log2(P_i), isJT);
        %[P_i2, EXITFLAG2, l2] = dinkelbach_algorithm(Pt, BW, w, R, R_min_JT_user, gamma, rho, P_fix, kappa, PCM, a_ib, c_ib, c1_b);

    end
  
%     P_i = 2.^q_test;
%     if(cvx_status=="Failed")
%         EXITFLAG = -8;
%     elseif( sum(A*P_i<=b) ~= length(b) )
%         EXITFLAG = -7;
%     else
%         EXITFLAG = 1;
%     end  
end

end