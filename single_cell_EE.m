function [EE, INUI_i, throughput_single_cell, p_tot_single_cell, user_i_data] = single_cell_EE(h_i, w, BW, P_i, ICI_i, rho, kappa, P_fix)
    N_users = length(h_i);
    
    
    % inter-NOMA-user interference (INUI)
    INUI_i = zeros(N_users,1);
    for k = 1:N_users
        for j = 1:(k-1)
            INUI_i(k) = INUI_i(k) + P_i(j).*h_i(k);
        end
    end
    
    % Throughput
    user_i_data = w*BW*log2(1 + (h_i.*P_i)./(ICI_i + INUI_i + w));
    throughput_single_cell = sum(w*BW*log2(1 + (h_i.*P_i)./(ICI_i + INUI_i + w)));
    
    %Power consumption
    p_tot_single_cell = 0;
    for i = 1:N_users
        p_tot_single_cell = p_tot_single_cell + P_i(i) + P_i(i).*rho + kappa*(N_users-i+1); 
    end
    p_tot_single_cell = p_tot_single_cell + P_fix;
    
    EE = throughput_single_cell/p_tot_single_cell;
end