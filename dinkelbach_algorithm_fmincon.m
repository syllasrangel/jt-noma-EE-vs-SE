function [P_i, EXITFLAG, l] = dinkelbach_algorithm_fmincon(Pt, BW, w, R, R_min_JT_user, gamma, rho, P_fix, kappa, PCM, a_ib, c_ib, c1_b, q_initial, isJT)
L_max=30; % iteration number
%isJT = true;

[N_users, N_BSs, ~] = size(gamma);
N_inner_users = N_users - 1;

A=[];
for bs = 1:N_BSs
   A = [A, [zeros(bs-1,N_users); ones(1, N_users); zeros(N_BSs-bs,N_users)]];
end
b = Pt*ones(N_BSs,1);


lambda = 0;


gamma_ib = NaN(N_users,N_BSs);
for bs = 1:N_BSs
    for i = 1:N_users
        gamma_ib(i,bs) = gamma(i,bs,bs);
    end
end

% Calculates the number of users per cluster
if(~isJT)
    J_b = N_inner_users*ones(N_BSs,1);
    J_b(1) = N_users;
else
    J_b = N_users*ones(N_BSs,1);
end

q_test = q_initial;


options = optimoptions('fmincon','Display','off','MaxFunctionEvaluations', 300000, 'MaxIterations', 40000, 'StepTolerance',1e-30, 'Algorithm','sqp');



fileID = fopen('dinkelbach_output.txt','w');
fprintf(fileID,'-----------------------\n');
if(isJT)
    fprintf(fileID,'JT-CoMP NOMA (fmincon)\n');
else
    fprintf(fileID,'Conventional NOMA (fmincon)\n');
end
fprintf(fileID,'a_ib:\n');
fprintf(fileID,'%f | %f \n',a_ib)
fprintf(fileID,'c_ib:\n');
fprintf(fileID,'%f | %f \n',c_ib)
fprintf(fileID,'c1_b:\n');
fprintf(fileID,'%f\n',c1_b)
fprintf(fileID,'-----------------------\n');

fprintf(fileID,'%-8s %-14s %-14s %-14s %-14s %-14s %-14s %-14s %-14s %-14s %-14s %-14s   difference in percent\n', 'iter', 'lambda (EE)', 'Data rate tot', 'P_11','P_21','P_31','P_12','P_22','P_32', 'fmincon_optval', 'calc_obj', 'optval_gap');

for l = 0:L_max
    
    obj_func = @(q)-EE_obj_function_global(w, BW, gamma, rho, P_fix, kappa, lambda, PCM, isJT, a_ib, c_ib, c1_b, q);
    
    nlcon = @(x)nonlcon_ICI_global(w, BW, gamma, R, R_min_JT_user, Pt, isJT, false, c1_b, x);   
    [q_test,fmincon_optval,EXITFLAG] = fmincon(obj_func,q_test,[],[],[],[],[],[],nlcon,options);
    
    [calc_obj, N, D] = EE_obj_function_global(w, BW, gamma, rho, P_fix, kappa, lambda, PCM, isJT, a_ib, c_ib, c1_b, q_test);
    %N-lambda*D
    %TODO: is it really abs?
    optval_gap = -fmincon_optval - calc_obj;
    [R_tot,~] = system_throughput(w, BW, gamma, isJT, 2.^q_test);
    calc_obj
    fprintf(fileID,'%-8d %-14.4f %-14.4f %-14.8f %-14.8f %-14.8f %-14.8f %-14.8f %-14.8f %-14.4f %-14.4f %-14.4f   %f%%\n', l, lambda, R_tot, 2.^q_test.', -fmincon_optval, calc_obj, optval_gap, -(optval_gap/fmincon_optval)*100);
%     if(abs(N-lambda*D) < 0.00001)
%         break;
%     end
    
    lambda = N/D;
    
end
  
fclose(fileID);
P_i = 2.^q_test;
% if(cvx_status=="Failed")
%     EXITFLAG = -8;
% elseif( sum(A*P_i<=b) ~= length(b) )
%     EXITFLAG = -7;
% else
%     EXITFLAG = 1;
% end  

end