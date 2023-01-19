function [P_i, EXITFLAG] = min_power_global_fmincon(Pt, BW, w, R, R_min_JT_user, gamma, rho, P_fix, kappa, PCM, isJT)
[N_users, N_BSs, ~] = size(gamma);
N_inner_users = N_users - 1;



%options = optimoptions('fmincon','Display','off','MaxFunctionEvaluations', 300000, 'MaxIterations', 40000, 'StepTolerance',1e-30);
options = optimoptions('fmincon','Display','off','MaxFunctionEvaluations', 300000, 'MaxIterations', 40000, 'StepTolerance',1e-30, 'Algorithm','sqp');

fun2 = @(x)sum(2.^x);

if (~isJT)
    L_max = 0;
else
    L_max = 50;
end
% p = 2^q => q = log2(p)
q_test = log2((Pt/N_users)*ones(N_BSs*(N_inner_users + isJT) + ~isJT,1));

c1_b = zeros(N_BSs,1);
c1_b(N_BSs) = 1;



% Calculates the number of users per cluster
if(~isJT)
    J_b = N_inner_users*ones(N_BSs,1);
    J_b(1) = N_users;
else
    J_b = N_users*ones(N_BSs,1);
end

gamma_ib = NaN(N_users,N_BSs);
for bs = 1:N_BSs
    for i = 1:N_users
        gamma_ib(i,bs) = gamma(i,bs,bs);
    end
end

for l = 0:L_max
    
    if(l~=0)
         %abs(system_throughput(w, BW, gamma, isJT, 2.^q_test)/system_power_consumption(2.^q_test, gamma, rho, P_fix, kappa, PCM, isJT) - system_throughput(w, BW, gamma, isJT, P_i)/system_power_consumption(P_i, gamma, rho, P_fix, kappa, PCM, isJT))
        if(abs(system_throughput(w, BW, gamma, isJT, 2.^q_test)/system_power_consumption(2.^q_test, gamma, rho, P_fix, kappa, PCM, isJT) - system_throughput(w, BW, gamma, isJT, P_i)/system_power_consumption(P_i, gamma, rho, P_fix, kappa, PCM, isJT)) < 0.0001)
            break;
        end
    end
    
    P_i = 2.^q_test;
    P_ib = Pvec2mat(gamma, isJT, P_i);

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

    if(isnan(c1_b))
        test = 1;
    end

    %c1_b
    nlcon = @(x)nonlcon_ICI_global(w, BW, gamma, R, R_min_JT_user, Pt, isJT, true, c1_b, x);   
    [q_test,~,EXITFLAG] = fmincon(fun2,q_test,[],[],[],[],[],[],nlcon,options);

    
end
P_i = 2.^q_test;


end