function [f, individual_rate] = system_throughput(w, BW, gamma, isJT, P_i)
    [N_users, N_BSs, ~] = size(gamma);
    N_inner_users = N_users - 1;
    
    individual_rate = NaN(N_users,N_BSs);
    
    % Reshapes the vector P_i   
    P_ib = Pvec2mat(gamma, isJT, P_i); 

    
    [ICI, INUI] = interference(gamma, isJT, P_i);
    f = 0;
    % Inner users data rate
    for bs = 1:N_BSs
        for i = 1:N_inner_users
            individual_rate(i,bs) = w*BW*log2(1 + (gamma(i,bs,bs)*P_ib(i,bs))./(ICI(i,bs) + INUI(i,bs) + w));
            f = f + individual_rate(i,bs);
        end
    end
    
    % Edge user data rate
    if(~isJT)
        individual_rate(N_users,1) = w*BW*log2(1 + (gamma(N_users,1,1)*P_ib(N_users,1))./(ICI(N_users,1) + INUI(N_users,1) + w));
        SINR = 10*log10((gamma(N_users,1,1)*P_ib(N_users,1))./(ICI(N_users,1) + INUI(N_users,1) + w));
        individual_rate(N_users,1);
    else
        sum_useful_term = 0;
        for bs_aux = 1:N_BSs
            sum_useful_term = sum_useful_term + P_ib(N_users,bs_aux).*gamma(N_users,bs_aux,bs_aux);
        end
        individual_rate(N_users,1) = w*BW*log2(1 + (sum_useful_term)./(ICI(N_users,1) + INUI(N_users,1) + w));
    end
    f = f + individual_rate(N_users,1);
end