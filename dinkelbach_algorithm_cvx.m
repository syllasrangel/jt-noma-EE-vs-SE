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



% fileID = fopen('dinkelbach_output.txt','w');
% fprintf(fileID,'-----------------------\n');
% if(isJT)
%     fprintf(fileID,'JT-CoMP NOMA\n');
% else
%     fprintf(fileID,'Conventional NOMA\n');
% end
% fprintf(fileID,'a_ib:\n');
% fprintf(fileID,'%f | %f \n',a_ib)
% fprintf(fileID,'c_ib:\n');
% fprintf(fileID,'%f | %f \n',c_ib)
% fprintf(fileID,'c1_b:\n');
% fprintf(fileID,'%f\n',c1_b)
% fprintf(fileID,'-----------------------\n');
% 
% fprintf(fileID,'%-8s %-14s %-14s %-14s %-14s %-14s %-14s %-14s %-14s %-14s %-14s %-14s   difference in percent\n', 'iter', 'lambda (EE)', 'Data rate tot', 'P_11','P_21','P_31','P_12','P_22','P_32', 'cvx_optval', 'calc_obj', 'cvx_optval_gap');

for l = 0:L_max

    cvx_begin quiet
        variable q_test(length_q);
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
            
            % === SIC constraint ===
            if(SIC_constraint)
                for bs = 1:N_BSs
                    for j = 1:(J_b(bs)-1)
                        for k = j+1:J_b(bs)
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
    
    % TODO: Save it to a output file in form of column with all the variable values for each iteration per line. 
    %[ttttt, N, D] = EE_obj_function_global(w, BW, gamma, rho, P_fix, kappa, lambda, PCM, isJT, a_ib, c_ib, c1_b, q_test);
    [calc_obj, N, D] = EE_obj_function_global_CVX(w, BW, gamma, rho, P_fix, kappa, lambda, PCM, isJT, a_ib, c_ib, c1_b, zeros(N_inner_users+1, N_BSs), zeros(N_inner_users+1, N_BSs), q_test);
    %calc_obj
    %cvx_optval
    %cvx_optval_gap = cvx_optval - calc_obj;
    %delta_opt_fraction = delta_opt/cvx_optval
    %TODO: is it really abs?
    %lambda
    
%     [R_tot,~] = system_throughput(w, BW, gamma, isJT, 2.^q_test);
    
%     fprintf(fileID,'%-8d %-14.4f %-14.4f %-14.8f %-14.8f %-14.8f %-14.8f %-14.8f %-14.8f %-14.4f %-14.4f %-14.4f   %f%%\n', l, lambda, R_tot, 2.^q_test.', cvx_optval, calc_obj, cvx_optval_gap, (cvx_optval_gap/cvx_optval)*100);
    if(N-lambda*D < 0.00001)
        break;
    end

    lambda = N/D;
    
    %lambda2 = N2/D2;
    
end

% fclose(fileID);

status = cvx_status;
P_i = 2.^q_test;

end