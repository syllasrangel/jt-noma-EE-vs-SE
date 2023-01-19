function [ICI, INUI] = interference_v2(gamma, isJT, P_ib)
    [N_users, N_BSs, ~] = size(gamma);
    N_inner_users = N_users - 1;
    
    % Calculates the number of users per cluster
    if(~isJT)
        J_b = N_inner_users*ones(N_BSs,1);
        % finds the bs with power different of NaN in to the edge user
        bs_edge_aux = find(~isnan(P_ib(N_users,:)));
        bs_edge = bs_edge_aux(1);
        J_b(bs_edge) = N_users;
    else
        J_b = N_users*ones(N_BSs,1);
    end
    
    % Reshapes the vector q
%     if(~isJT)% NaN in q(JT_user,: ~1) for conventional NOMA
%         insert = @(a, x, n)cat(1,  x(1:n), a, x(n+1:end));
%         for bs1 = 2:N_BSs
%             %q = insert(-inf, q, bs1*(N_inner_users + 1) - 1);
%             P_i = insert(NaN, P_i, bs1*(N_inner_users + 1) - 1);
%         end
%     end
%     P_ib = reshape(P_i,[N_users,N_BSs]);
    
    % ICI calculation for both conventional and JT NOMA (depending on 'isJT')
    ICI = zeros(N_inner_users+1, N_BSs);
    if(~isJT) % NaN in ICI(JT_user,: ~1) for conventional NOMA
        ICI(N_users,:) = NaN;
        ICI(N_users, bs_edge) = 0;
    end
    for bs1 = 1:N_BSs
        for k = 1:J_b(bs1)
            for bs2 = 1:N_BSs
                if(bs1~=bs2)
                    %for j = 1:(J_b(bs2) - isJT) % When JT it does not consider interference from the cell-edge user
                    for j = 1:(J_b(bs2) - ( isJT && k == J_b(bs2) ) ) % When JT it does not consider interference from the cell-edge user on itself
                        ICI(k,bs1) = ICI(k,bs1) + P_ib(j,bs2).*gamma(k,bs1,bs2);
                    end
                end
            end
        end
    end
    
    % inter-NOMA-user interference (INUI)
    INUI = zeros(N_inner_users+1, N_BSs);
    if(~isJT) % NaN in INUI(JT_user,: ~1) for conventional NOMA
        INUI(N_users,:) = NaN;
        INUI(N_users, bs_edge) = 0;
    end
    for bs1 = 1:N_BSs
        for k = 1:J_b(bs1)
            for j = 1:(k-1)
                INUI(k,bs1) = INUI(k,bs1) + P_ib(j,bs1).*gamma(k,bs1,bs1);
            end
        end
    end
    
    
end