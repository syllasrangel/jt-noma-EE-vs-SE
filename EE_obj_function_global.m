function [f, N, D] = EE_obj_function_global(w, BW, gamma, rho, P_fix, kappa, lambda, PCM, isJT, a_ib, c_ib, c1_b, q)
    [N_users, N_BSs, ~] = size(gamma);
    N_inner_users = N_users - 1;    

%     N = system_throughput(w, BW, gamma, isJT, P_i);

    % Reshapes the vector q
    q_ib = Pvec2mat(gamma, isJT, q);    
    
    % Calculates numerator
    P_i = 2.^q;
    [ICI, INUI] = interference(gamma, isJT, P_i);
    
    
    N = 0;
    for bs = 1:N_BSs
        for i = 1:N_inner_users
            N = N + a_ib(i,bs)*w*BW*(log2(gamma(i,bs,bs)) + q_ib(i,bs)) + c_ib(i,bs)*w*BW - a_ib(i,bs)*w*BW*log2(ICI(i,bs) + INUI(i,bs) + w);
        end
    end
    
    % Edge user data rate
    if(~isJT)
         N = N + a_ib(N_users,1)*w*BW*(log2(gamma(N_users,1,1)) + q_ib(N_users,1)) + c_ib(N_users,1)*w*BW - a_ib(N_users,1)*w*BW*log2(ICI(N_users,1) + INUI(N_users,1) + w);
    else
        sum_aux= 0;
        for bs_aux = 1:N_BSs
            if(c1_b(bs_aux)~=0) % Workaround for CVX error ({real affine} + {invalid})
                sum_aux = sum_aux + c1_b(bs_aux)*(q_ib(N_users,bs_aux) + log2(gamma(N_users,bs_aux,bs_aux)/c1_b(bs_aux)));
            end
            %sum_aux = sum_aux + c1_b(bs_aux)*q_ib(N_users,bs_aux) + log2((gamma(N_users,bs_aux,bs_aux)/c1_b(bs_aux)).^c1_b(bs_aux));
        end 
        N = N + a_ib(N_users,1)*w*BW*sum_aux + c_ib(N_users,1)*w*BW - a_ib(N_users,1)*w*BW*log2(ICI(N_users,1) + INUI(N_users,1) + w);
    end
    
    % Calculates denominator
    D = system_power_consumption(P_i, gamma, rho, P_fix, kappa, PCM, isJT);

    f = N - lambda*D;

end