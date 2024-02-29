function [P_i, status, l] = dinkelbach_algorithm_local_cvx(Pt, BW, w, R_min_i, h_i, rho, P_fix, kappa, a_i, c_i, ICI_i)
L_max=30; % Max number of iterations

N_users = length(h_i);


lambda = 0;


for l = 0:L_max

    cvx_begin quiet
        variable q(N_users);
        expression INUI_i(N_users)
        maximize( EE_obj_function_local_cvx(w, BW, h_i, rho, kappa, P_fix, a_i, c_i, ICI_i, INUI_i, lambda, q));
        subject to
            for k = 1:N_users
                INUI_i(k)=0;
                for j = 1:(k-1)
                    INUI_i(k) = INUI_i(k) + (2.^q(j)).*h_i(k);
                end
            end

            % === Rate requirement inner users ===
            for j = 1:N_users
                gamma_min = 2^(R_min_i(j)/(w*BW)) - 1;
                -(q(j) - log(ICI_i(j) + INUI_i(j) + w)./log(2) + log(h_i(j)/gamma_min)./log(2)) <= 0;
            end

            % === BS max power constraint ===
            tot_power = 0;
            for j = 1:N_users 
                tot_power = tot_power + 2.^(q(j));
            end
            tot_power - Pt <= 0;
            
    cvx_end
    
    [~, N, D] = EE_obj_function_local_cvx(w, BW, h_i, rho, kappa, P_fix, a_i, c_i, ICI_i, zeros(N_users,1), lambda, q);
    
    if(N-lambda*D < 0.00001 || status2exitflag(cvx_status)<0)
        break;
    end

    lambda = N/D;
end


status = cvx_status;
P_i = 2.^q;

end