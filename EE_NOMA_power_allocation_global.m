function [P_i, EXITFLAG, l] = EE_NOMA_power_allocation_global(Pt, P_tol, BW, w, R, R_min_JT_user, gamma, rho, P_fix, kappa, PCM, isJT)
% TODO: Ver se é melhor fazer tanto JT quanto convencional aqui ou usar funções separadas <---------------------------------------------------
L_max=30; % iteration number

%========= Initialization ==========
%options = optimoptions('fmincon','Display','off','MaxFunctionEvaluations', 300000, 'MaxIterations', 40000, 'StepTolerance',1e-30, 'UseParallel', true);%, 'Algorithm','sqp');
%options = optimoptions('fmincon','Display','off','MaxFunctionEvaluations', 3000000, 'MaxIterations', 400000, 'StepTolerance',1e-30, 'ConstraintTolerance', 1e-30, 'UseParallel', true);
%options = optimoptions('fmincon','Display','off','MaxFunctionEvaluations', 300000, 'MaxIterations', 40000, 'StepTolerance',1e-30, 'ConstraintTolerance', 1e-30);
%options = optimoptions('fmincon','Display','off','MaxFunctionEvaluations', 300000, 'MaxIterations', 40000, 'StepTolerance',1e-30, 'Algorithm','sqp');
%options = optimoptions('fmincon','Display','off','MaxFunctionEvaluations', 300000, 'MaxIterations', 40000, 'StepTolerance',1e-30);


[N_users, N_BSs, ~] = size(gamma);
N_inner_users = N_users - 1;


A=[];
for bs = 1:N_BSs
   aux = (bs ~= 1 && ~isJT);
   A = [A, [zeros(bs-1,N_users - aux); ones(1, N_users - aux); zeros(N_BSs-bs,N_users - aux)]];
end
b = Pt*ones(N_BSs,1);


% Initial power minimization problem
% p = 2^q => q = log2(p)
% x0 = log2((Pt/N_users)*ones(N_BSs*(N_inner_users + isJT) + ~isJT,1));
% if(isJT)
% x0 = [2.2400;
%     2.2402;
%     3.2351;
%     2.7094;
%     2.7577;
%     2.7361];
% else
%     x0 = [1.7655;
%     1.7698;
%     3.8280;
%     1.5279;
%     1.5345];
% end
% q0_ib = Pvec2mat(gamma, isJT, x0);
% c1_b = NaN(N_BSs,1);
% if(isJT)
%     sum_aux = 0;
%     for bs = 1:N_BSs
%         sum_aux = sum_aux + 2.^q0_ib(N_users,bs).*gamma(N_users,bs,bs);
%     end
%     for bs = 1:N_BSs
%         c1_b(bs) = 2.^q0_ib(N_users,bs).*gamma(N_users,bs,bs)./sum_aux;
%     end
% end
% nlcon = @(x)nonlcon_ICI_global(w, BW, gamma, R, R_min_JT_user, Pt, P_tol, isJT, true, c1_b, x);
 fun2 = @(x)sum(2.^x);
% funZero = @(x)0;
% [P_initial_zero,~,EXITFLAG1_zero] = fmincon(funZero,x0,[],[],[],[],[],[],nlcon,options);
% [P_initial,~,EXITFLAG1] = fmincon(fun2,P_initial_zero,[],[],[],[],[],[],nlcon,options);
% 
% q_i = P_initial;

lambda = 0;

length_q = N_BSs*(N_inner_users + isJT) + ~isJT;

% cvx_begin quiet
%     variable q_test(length_q);
%     expression INUI(N_inner_users+1, N_BSs)
%     expression ICI(N_inner_users+1, N_BSs)
%     minimize( fun2(q_test) );
%     subject to
%         % Calculates the number of users per cluster
%         if(~isJT)
%             J_b = N_inner_users*ones(N_BSs,1);
%             J_b(1) = N_users;
%         else
%             J_b = N_users*ones(N_BSs,1);
%         end
% 
% %         % Reshapes the vector q   
% %         q_ib = Pvec2mat(gamma, isJT, q_test);
% 
%         % ICI and inter-NOMA-user interference (INUI) calculation
%         [ICI, INUI] = interference_CVX(gamma, isJT, 2.^q_test, ICI, INUI);
% 
% %         % === SIC constraint ===
% %         for bs = 1:N_BSs
% %             for k = 1:J_b(bs)-1
% %                 for ii = k+1:J_b(bs)
% %                     sum_INUI_power = 0;
% %                     for aa = 1:ii-1
% %                         sum_INUI_power = sum_INUI_power + 2.^q_test(two_dim_2_one_dim(aa,bs,N_users, isJT)).*gamma(k,bs,bs);
% %                     end
% %                     if(isJT && ii == J_b(bs))
% %                         sum_useful_power = 0;
% %                         for bs_aux = 1:N_BSs
% %                             sum_useful_power = sum_useful_power + 2.^q_test(two_dim_2_one_dim(ii,bs_aux,N_users, isJT)).*gamma(k,bs_aux,bs_aux);
% %                         end
% %                          -(sum_useful_power - sum_INUI_power - ICI(k,bs) - P_tol) <= 0;
% %                     else
% %                         %sum_useful_power = 2.^q_test(two_dim_2_one_dim(ii,bs,N_users, isJT)).*gamma(k,bs,bs);
% %                         %New formulation
% %                         -q_test(two_dim_2_one_dim(ii,bs,N_users, isJT)) - log(gamma(k,bs,bs))./log(2) + log(sum_INUI_power + ICI(k,bs) + P_tol)./log(2) <= 0
% %                     end
% %                     %idx_counter = idx_counter+1;
% %                     %-(sum_useful_power - sum_INUI_power - ICI(k,bs) - P_tol) <= 0;
% %                 end
% %             end
% %         end
% 
%         % === Rate requirement inner users ===
%         for bs = 1:N_BSs
%             for j = 1:N_inner_users
%                 gamma_min = 2^(R(j,bs)/(w*BW)) - 1;
%                 %c(idx_counter) = -(q_ib(j,bs) - log2(ICI(j,bs) + INUI(j,bs) + w) + log2(gamma(j,bs,bs)/gamma_min));
%                 -(q_test(two_dim_2_one_dim(j, bs, N_users, isJT)) - log(ICI(j,bs) + INUI(j,bs) + w)./log(2) + log(gamma(j,bs,bs)/gamma_min)./log(2)) <= 0;
%             end
%         end
% 
%         % === Rate requirement cell-edge user ===
%         gamma_min = 2^(R_min_JT_user/(w*BW)) - 1;
%         if(~isJT)
%             %c(idx_counter) = -(q_ib(J_b(1),1) - log2(ICI(J_b(1),1) + INUI(J_b(1),1) + w) + log2(gamma(J_b(1),1,1)/gamma_min));
%             -(q_test(two_dim_2_one_dim(J_b(1),1,N_users, isJT)) - log(ICI(J_b(1),1) + INUI(J_b(1),1) + w)./log(2) + log(gamma(J_b(1),1,1)/gamma_min)./log(2))<=0;
%         else
%             sum_useful_term = 0;
%             for bs_aux = 1:N_BSs % Relaxed constraint for JT case.
%                 sum_useful_term = sum_useful_term + q_test(two_dim_2_one_dim(J_b(bs_aux), bs_aux, N_users, isJT)).*c1_b(bs_aux) + c1_b(bs_aux).*log2(gamma(J_b(bs_aux),bs_aux,bs_aux)/c1_b(bs_aux));
%             end 
%             %c(idx_counter) = -(sum_useful_term - log2(ICI(J_b(1),1) + INUI(J_b(1),1) + w) - log2(gamma_min));
%             -(sum_useful_term - log(ICI(J_b(1),1) + INUI(J_b(1),1) + w)./log(2) - log(gamma_min)./log(2))<=0;
%         end
%         
% cvx_end


[P_i, ~] = min_power_global(Pt, BW, w, R, R_min_JT_user, gamma, rho, P_fix, kappa, PCM, isJT);

q_test = log2(P_i);


% Calculates the number of users per cluster
if(~isJT)
    J_b = N_inner_users*ones(N_BSs,1);
    J_b(1) = N_users;
else
    J_b = N_users*ones(N_BSs,1);
end

%if(EXITFLAG1 == -2 || (EXITFLAG1 == 1 && sum(A*(2.^P_initial)<=b) ~= length(b)))
if( sum(A*(2.^q_test)<=b) ~= length(b) )
    EXITFLAG = -6;
    l=0;
    P_i = 2.^q_test;
else
    for l = 0:L_max

        %---------------------------========================    
        %[~, N, D] = EE_obj_function_ICI_global(w, BW, gamma1_i, gamma2_i, gamma_2_1, gamma_1_2, rho, P_fix, kappa, lambda, PCM, isJT, P_i);
%         [~, N, D] = EE_obj_function_global(w, BW, gamma, rho, P_fix, kappa, lambda, PCM, isJT, a_ib, c_ib, q_i);
%         nlcon = @(x)nonlcon_ICI_global(w, BW, gamma, R, R_min_JT_user, Pt, P_tol, isJT, false, c_b, x);
        %---------------------------========================
%         if(l>10)
%             test=1;
%         end
        
        if(l~=0)
            %if(abs(N-lambda*D)<0.0001)
            %if(abs(system_throughput(w, BW, gamma, isJT, 2.^q_i)/system_power_consumption(2.^q_i, gamma, rho, P_fix, kappa, PCM, isJT) - system_throughput(w, BW, gamma, isJT, P_i)/system_power_consumption(P_i, gamma, rho, P_fix, kappa, PCM, isJT)) < 0.0001)
            %abs(system_throughput(w, BW, gamma, isJT, 2.^q_test)/system_power_consumption(2.^q_test, gamma, rho, P_fix, kappa, PCM, isJT) - system_throughput(w, BW, gamma, isJT, P_i)/system_power_consumption(P_i, gamma, rho, P_fix, kappa, PCM, isJT))
            if(abs(system_throughput(w, BW, gamma, isJT, 2.^q_test)/system_power_consumption(2.^q_test, gamma, rho, P_fix, kappa, PCM, isJT) - system_throughput(w, BW, gamma, isJT, P_i)/system_power_consumption(P_i, gamma, rho, P_fix, kappa, PCM, isJT)) < 0.0001)
                break;
            end
            
%             % Used to save time in outage cases. TODO
%             if(EXITFLAG1 ~= 1 && EXITFLAG ~= 1 && l > 15)
%                 break;
%             end
        end
        P_i = 2.^q_test;
        
        %SINR for non-CoMP case
        [ICI, INUI] = interference(gamma, isJT, P_i);
        P_ib = Pvec2mat(gamma, isJT, P_i);
        
        gamma_ib = NaN(N_users,N_BSs);
        for bs = 1:N_BSs
            for i = 1:N_users
                gamma_ib(i,bs) = gamma(i,bs,bs);
            end
        end
        SINR_ib = P_ib.*gamma_ib./(ICI + INUI + w);
    
        a_ib = SINR_ib./(1+SINR_ib);
        c_ib = log2(1+SINR_ib) - SINR_ib.*log2(SINR_ib)./(1+SINR_ib);
        if(isnan(a_ib(1,1)) || isnan(c_ib(1,1))) % TODO: remove if
            test = 1
        end
        
        c1_b = NaN(N_BSs,1);
        if(isJT)
            sum_aux = 0;
            for bs = 1:N_BSs
                sum_aux = sum_aux + P_ib(N_users,bs).*gamma_ib(N_users,bs);
            end
            for bs = 1:N_BSs
                c1_b(bs) = P_ib(N_users,bs).*gamma_ib(N_users,bs)./sum_aux;
            end
        end
        
        [~, N, D] = EE_obj_function_global(w, BW, gamma, rho, P_fix, kappa, lambda, PCM, isJT, a_ib, c_ib, c1_b, q_test);
        
        lambda = N/D;
        
        if(isnan(lambda))
            test=1;
        end
        
%         nlcon = @(x)nonlcon_ICI_global(w, BW, gamma, R, R_min_JT_user, Pt, P_tol, isJT, false, c1_b, x);
%         fun = @(x)-EE_obj_function_global(w, BW, gamma, rho, P_fix, kappa, lambda, PCM, isJT, a_ib, c_ib, c1_b, x);
% 
% 
%         [q_i,~,EXITFLAG] = fmincon(fun,q_i,[],[],[],[],[],[],nlcon,options);

        cvx_begin quiet
            variable q_test(length_q);
            expression INUI(N_inner_users+1, N_BSs)
            expression ICI(N_inner_users+1, N_BSs)
            minimize( -EE_obj_function_global_CVX(w, BW, gamma, rho, P_fix, kappa, lambda, PCM, isJT, a_ib, c_ib, c1_b, INUI, ICI, q_test) );
            subject to
                % Calculates the number of users per cluster
                if(~isJT)
                    J_b = N_inner_users*ones(N_BSs,1);
                    J_b(1) = N_users;
                else
                    J_b = N_users*ones(N_BSs,1);
                end

        %         % Reshapes the vector q   
        %         q_ib = Pvec2mat(gamma, isJT, q_test);

                % ICI and inter-NOMA-user interference (INUI) calculation
                [ICI, INUI] = interference_CVX(gamma, isJT, 2.^q_test, ICI, INUI);


%                 % === SIC constraint ===
%                 for bs = 1:N_BSs
%                     for k = 1:J_b(bs)-1
%                         for ii = k+1:J_b(bs)
%                             sum_INUI_power = 0;
%                             for aa = 1:ii-1
%                                 sum_INUI_power = sum_INUI_power + 2.^q_test(two_dim_2_one_dim(aa,bs,N_users, isJT)).*gamma(k,bs,bs);
%                             end
%                             if(isJT && ii == J_b(bs))
%                                 sum_useful_power = 0;
%                                 for bs_aux = 1:N_BSs
%                                     sum_useful_power = sum_useful_power + 2.^q_test(two_dim_2_one_dim(ii,bs_aux,N_users, isJT)).*gamma(k,bs_aux,bs_aux);
%                                 end
%                                  -(sum_useful_power - sum_INUI_power - ICI(k,bs) - P_tol) <= 0;
%                             else
%                                 %sum_useful_power = 2.^q_test(two_dim_2_one_dim(ii,bs,N_users, isJT)).*gamma(k,bs,bs);
%                                 %New formulation
%                                 -q_test(two_dim_2_one_dim(ii,bs,N_users, isJT)) - log(gamma(k,bs,bs))./log(2) + log(sum_INUI_power + ICI(k,bs) + P_tol)./log(2) <= 0
%                             end
%                             %idx_counter = idx_counter+1;
%                             %-(sum_useful_power - sum_INUI_power - ICI(k,bs) - P_tol) <= 0;
%                         end
%                     end
%                 end

                % === Rate requirement inner users ===
                for bs = 1:N_BSs
                    for j = 1:N_inner_users
                        gamma_min = 2^(R(j,bs)/(w*BW)) - 1;
                        %c(idx_counter) = -(q_ib(j,bs) - log2(ICI(j,bs) + INUI(j,bs) + w) + log2(gamma(j,bs,bs)/gamma_min));
                        -(q_test(two_dim_2_one_dim(j, bs, N_users, isJT)) - log(ICI(j,bs) + INUI(j,bs) + w)./log(2) + log(gamma(j,bs,bs)/gamma_min)./log(2)) <= 0;
                    end
                end

                % === Rate requirement cell-edge user ===
                gamma_min = 2^(R_min_JT_user/(w*BW)) - 1;
                if(~isJT)
                    %c(idx_counter) = -(q_ib(J_b(1),1) - log2(ICI(J_b(1),1) + INUI(J_b(1),1) + w) + log2(gamma(J_b(1),1,1)/gamma_min));
                    -(q_test(two_dim_2_one_dim(J_b(1),1,N_users, isJT)) - log(ICI(J_b(1),1) + INUI(J_b(1),1) + w)./log(2) + log(gamma(J_b(1),1,1)/gamma_min)./log(2))<=0;
                else
                    sum_useful_term = 0;
                    for bs_aux = 1:N_BSs % Relaxed constraint for JT case.
                        sum_useful_term = sum_useful_term + q_test(two_dim_2_one_dim(J_b(bs_aux), bs_aux, N_users, isJT)).*c1_b(bs_aux) + c1_b(bs_aux).*log(gamma(J_b(bs_aux),bs_aux,bs_aux)/(log(2)*c1_b(bs_aux)));
                    end 
                    %c(idx_counter) = -(sum_useful_term - log2(ICI(J_b(1),1) + INUI(J_b(1),1) + w) - log2(gamma_min));
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
        cvx_end

    end
    
    P_i = 2.^q_test;
    %2.^q_i


    if( sum(A*P_i<=b) ~= length(b) )
        EXITFLAG = -7;
    else
        EXITFLAG = 1;
    end  
end



% if(EXITFLAG ~= 1 && sum(A*P_i<=b) ~= length(b) )
%     EXITFLAG = -7;
% end

end