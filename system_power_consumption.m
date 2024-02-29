function Pt = system_power_consumption(P_i, gamma, rho, P_fix, kappa, isJT, is_cvx)
    [N_users, N_BSs, ~] = size(gamma);
    N_inner_users = N_users - 1;
    
    force_PCM_conv = true; % forces the PCM conv in the cvx allocation

       
    %P_ib = Pvec2mat(gamma, isJT, P_i);
    P_ib = Pvec2mat(gamma, length(P_i) == N_BSs*N_users, P_i);
    
    Pt = 0;
    for bs = 1:N_BSs
        Pt = Pt + P_fix + sum(user_power_consumption_2(P_ib(:,bs), gamma, rho, kappa, is_cvx));
%         if(is_cvx)
%             Pt = Pt + xb(bs)*((N_users-1)*kappa); %new: adds an additional decoding for each non-CoMP user case x = 1, i.e., power for CoMP user > 0            
%         end
    end
    
%     if(is_cvx)
%         Pt = Pt + kappa; % new: power for the edge user decode its own signal
%     end
    
    % Workaround: for PCM proposed in CVX. 
    % If isJT is true and force_PCM_conv false, it will calculate the SIC
    % overhead for edge user in all BS even if the allocated power is zero
    % in some BS.
% NEW: Uncomment this, case remove xb
    if(is_cvx && ~force_PCM_conv && isJT)
        Pt = Pt + N_BSs*(N_users -1)*kappa + kappa;
    elseif(is_cvx && (~isJT || force_PCM_conv)) % if is CVX and conv. then only one BS transmit for the edge user.
        Pt = Pt + N_users*kappa;
    end
    
end