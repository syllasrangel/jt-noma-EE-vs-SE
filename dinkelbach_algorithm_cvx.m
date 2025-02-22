function [P_i, status, l] = dinkelbach_algorithm_cvx(Pt, BW, w, R, R_min_JT_user, gamma, rho, P_fix, kappa, PCM, a_ib, c_ib, c1_b, isJT, SIC_constraint)
%SIC_constraint = true;
L_max=30; % Max number of iterations

[N_users, N_BSs, ~] = size(gamma);
N_inner_users = N_users - 1;

% Calculates the number of users per cluster
if(~isJT)
    J_b = N_inner_users*ones(N_BSs,1);
    J_b(1) = N_users;
else
    J_b = N_users*ones(N_BSs,1);
end


lambda = 0;

length_q = N_BSs*(N_inner_users + isJT) + ~isJT;


for l = 0:L_max

    cvx_begin quiet
        variable q_test(length_q);
        %variable x(N_BSs) % New variable
        %variable z(N_BSs) % New variable
        expression INUI(N_inner_users+1, N_BSs)
        expression ICI(N_inner_users+1, N_BSs)
        expression INUI_jk(N_inner_users+1, N_BSs, N_inner_users+1)
        maximize( EE_obj_function_global_CVX(w, BW, gamma, rho, P_fix, kappa, lambda, PCM, isJT, a_ib, c_ib, c1_b, INUI, ICI, q_test) );
        subject to
            % ICI and inter-NOMA-user interference (INUI) calculation
            [ICI, INUI, INUI_jk] = interference_CVX(gamma, isJT, 2.^q_test, ICI, INUI, INUI_jk);

            % === Rate requirement inner users ===
            for bs = 1:N_BSs
                for j = 1:N_inner_users
                    gamma_min = 2^(R(j,bs)/(w*BW)) - 1;
                    -(q_test(two_dim_2_one_dim(j, bs, N_users, isJT)) - log(ICI(j,bs) + INUI(j,bs) + w)./log(2) + log(gamma(j,bs,bs)/gamma_min)./log(2)) <= 0;
                end
            end

            % === Rate requirement cell-edge user ===
            gamma_min = 2^(R_min_JT_user/(w*BW)) - 1;               
            if(~isJT)
                -(q_test(two_dim_2_one_dim(J_b(1),1,N_users, isJT)) - log(ICI(J_b(1),1) + INUI(J_b(1),1) + w)./log(2) + log(gamma(J_b(1),1,1)/gamma_min)./log(2))<=0;
            else
                sum_useful_term = 0;
                for bs_aux = 1:N_BSs % Relaxed constraint for JT case.
                    sum_useful_term = sum_useful_term + q_test(two_dim_2_one_dim(J_b(bs_aux), bs_aux, N_users, isJT)).*c1_b(bs_aux) + log((gamma(J_b(bs_aux),bs_aux,bs_aux)./c1_b(bs_aux)).^c1_b(bs_aux))./log(2);
                end 
                -(sum_useful_term - log(ICI(J_b(1),1) + INUI(J_b(1),1) + w)./log(2) - log(gamma_min)./log(2))<=0;
            end

            % === BS max power constraint ===
            for bs = 1:N_BSs
                tot_power = 0;
                for j = 1:J_b(bs) 
                    tot_power = tot_power + 2.^(q_test(two_dim_2_one_dim(j, bs, N_users, isJT)));
                end
                tot_power - Pt <= 0;
            end
%             % === New constraints ===
%             if(isJT)
%                 x>=0;
%                 x<=1;
%                 z>=0;
%                 z<=Pt*x;
%                 for bs = 1:N_BSs
%                     -Pt*(1-x(bs)) <= z(bs) - 2.^q_test(two_dim_2_one_dim(J_b(1), bs, N_users, isJT))
%                     z(bs) <= 2.^q_test(two_dim_2_one_dim(J_b(1), bs, N_users, isJT))
%                 end
%                 CoMP_power = 0;
%                 for bs = 1:N_BSs
%                     CoMP_power = CoMP_power + 2.^q_test(two_dim_2_one_dim(J_b(1), bs, N_users, isJT));
%                 end
%                 sum(z) = CoMP_power;
%             end
            
            % === SIC constraint ===
            if(SIC_constraint)
                for bs = 1:N_BSs
                    %for j = 1:(J_b(bs)-1)
                        %for k = j+1:J_b(bs)
                    for j = 1:(J_b(bs)-1-isJT) % Not testing conditions to decode edge user`s signal in any user (it allows edge user to be served by only one BS)
                        for k = j+1:J_b(bs)-isJT % Not testing conditions to decode edge user`s signal in any user (it allows edge user to be served by only one BS)
                            if((isJT && k==J_b(bs)) || (~isJT && k==J_b(bs) && bs==1))
                                gamma_min = 2^(R_min_JT_user/(w*BW)) - 1;
                            else
                                gamma_min = 2^(R(k,bs)/(w*BW)) - 1;
                            end
                            -(q_test(two_dim_2_one_dim(k, bs, N_users, isJT)) - log(ICI(j,bs) + INUI_jk(j,bs,k) + w)./log(2) + log(gamma(j,bs,bs)/gamma_min)./log(2)) <= 0;
                        end
                    end
                end
            end
            
    cvx_end
    
    
    [calc_obj, N, D] = EE_obj_function_global_CVX(w, BW, gamma, rho, P_fix, kappa, lambda, PCM, isJT, a_ib, c_ib, c1_b, zeros(N_inner_users+1, N_BSs), zeros(N_inner_users+1, N_BSs), q_test);
    %calc_obj
    %cvx_optval
    %cvx_optval_gap = cvx_optval - calc_obj;
    %delta_opt_fraction = delta_opt/cvx_optval
    
    if(N-lambda*D < 0.00001)
        break;
    end

    lambda = N/D;
    
    
end


status = cvx_status;
P_i = 2.^q_test;

end