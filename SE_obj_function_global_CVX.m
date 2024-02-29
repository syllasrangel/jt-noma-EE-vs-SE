function N = SE_obj_function_global_CVX(w, BW, gamma, isJT, a_ib, c_ib, c1_b, INUI, ICI, q)
    [N_users, N_BSs, ~] = size(gamma);
    N_inner_users = N_users - 1;    


%     N = system_throughput(w, BW, gamma, isJT, P_i);


    two2one = @(j,bs)two_dim_2_one_dim(j, bs, N_users, isJT);   
    
    for bs=1:N_BSs
        for j = 1:N_inner_users+1
            ICI(j,bs)=0;
            INUI(j,bs)=0;
        end
    end
    
    % Calculates numerator
    P_i = 2.^q;
    [ICI, INUI] = interference_CVX(gamma, isJT, P_i, INUI, ICI);
    
    
    N = 0;
    for bs = 1:N_BSs
        for i = 1:N_inner_users
            N = N + a_ib(i,bs)*w*BW*(log(gamma(i,bs,bs))./log(2) + q(two2one(i,bs)) ) + c_ib(i,bs)*w*BW - a_ib(i,bs)*w*BW*log(ICI(i,bs) + INUI(i,bs) + w)./log(2);
        end
    end
    
    % Edge user data rate
    if(~isJT)
         N = N + a_ib(N_users,1)*w*BW*( log(gamma(N_users,1,1))./log(2) + q(two2one(N_users,1)) ) + c_ib(N_users,1)*w*BW - a_ib(N_users,1)*w*BW*log(ICI(N_users,1) + INUI(N_users,1) + w)./log(2);
    else
        sum_aux= 0;
        for bs_aux = 1:N_BSs
            if(c1_b(bs_aux)~=0) % Workaround for CVX error ({real affine} + {invalid})
                sum_aux = sum_aux + c1_b(bs_aux)*( q(two2one(N_users,bs_aux)) + log(gamma(N_users,bs_aux,bs_aux)/c1_b(bs_aux))./log(2) );
            end
        end 
        try
            a_ib(N_users,1)*w*BW*sum_aux + c_ib(N_users,1)*w*BW - a_ib(N_users,1)*w*BW*log(ICI(N_users,1) + INUI(N_users,1) + w)./log(2);
        catch
            error = true
        end
        N = N + a_ib(N_users,1)*w*BW*sum_aux + c_ib(N_users,1)*w*BW - a_ib(N_users,1)*w*BW*log(ICI(N_users,1) + INUI(N_users,1) + w)./log(2);
    end


end