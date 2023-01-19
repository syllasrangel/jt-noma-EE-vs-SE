function [P_i_tot, P_i_tot_user_side] = user_power_consumption_2(P_i, gamma, rho, kappa, is_cvx)
    [N_users, ~, ~] = size(gamma);

   
    P_i_tot = P_i + P_i*rho;
    for i=1:N_users
        % Workaround: if is_cvx assume that the edge-user user is not
        % served. The adittional power for that user is added in the
        % system_power_consumption
        if(is_cvx || (length(P_i)<N_users || P_i(N_users)==0)) % if there is no power allocated to the edge user
            P_i_tot(i) = P_i_tot(i) + kappa*(N_users-i);
        else
            P_i_tot(i) = P_i_tot(i) + kappa*(N_users-i+1);
        end
    end
    
    
    % User-side power
    P_i_tot_user_side = P_i*rho;
    %P_i_tot_user_side = P_i*0;
    for i=1:N_users
        % Workaround: if is_cvx assume that the edge-user user is not
        % served. The adittional power for that user is added in the
        % system_power_consumption
        if(is_cvx || (length(P_i)<N_users || P_i(N_users)==0)) % if there is no power allocated to the edge user
            P_i_tot_user_side(i) = P_i_tot_user_side(i) + kappa*(N_users-i);
        else
            P_i_tot_user_side(i) = P_i_tot_user_side(i) + kappa*(N_users-i+1);
        end
    end
end
