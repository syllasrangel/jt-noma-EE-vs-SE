function P_ib = Pvec2mat(gamma, isJT, P_i)
	[N_users, N_BSs, ~] = size(gamma);
    N_inner_users = N_users - 1;
    
    % Reshapes the vector P_i
    sizeP = size(P_i);
    if(sizeP(1)==1)
        P_i_aux = P_i.';
    else
        P_i_aux = P_i;
    end
    if(~isJT)% NaN in q(JT_user,: ~1) for conventional NOMA
        insert = @(a, x, n)cat(1,  x(1:n), a, x(n+1:end));
        for bs = 2:N_BSs
            %P_i_aux = insert(NaN, P_i_aux, bs*(N_inner_users + 1) - 1);
            P_i_aux = insert(0, P_i_aux, bs*(N_inner_users + 1) - 1);
        end
    end
    P_ib = reshape(P_i_aux,[N_users,N_BSs]);
end