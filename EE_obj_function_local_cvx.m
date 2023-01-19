function [f, N, D] = EE_obj_function_local_cvx(w, BW, h_i, rho, kappa, P_fix, a_i, c_i, ICI_i, INUI_i, lambda, q)
                                                
    N_users = length(h_i);


    % Calculates numerator
    P_i = 2.^q;
    
    %Power consumption
    D = 0;
    for i = 1:N_users
        D = D + P_i(i) + P_i(i).*rho + kappa*(N_users-i+1); 
    end
    D = D + P_fix;
    
    % inter-NOMA-user interference (INUI)
    for k = 1:N_users
        INUI_i(k)=0;
        for j = 1:(k-1)
            INUI_i(k) = INUI_i(k) + P_i(j).*h_i(k);
        end
    end    
    
    N = 0;
    for i = 1:N_users
        N = N + a_i(i)*w*BW*(log(h_i(i))./log(2) + q(i) ) + c_i(i)*w*BW - a_i(i)*w*BW*log(ICI_i(i) + INUI_i(i) + w)./log(2);
    end
    
    f = N - lambda*D;

end