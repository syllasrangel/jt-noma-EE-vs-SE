function [P_i, EXITFLAG, l, EE_total, P_BS1, P_BS2] = EE_NOMA_power_allocation_local_cvx(Pt, BW, w, R, R_min_JT_user, gamma, rho, P_fix, kappa, PCM)
L_max=50; % Max number of iterations
isJT = false;

[N_users, N_BSs, ~] = size(gamma);


gamma_ib = NaN(N_users,N_BSs);
for bs = 1:N_BSs
    for i = 1:N_users
        gamma_ib(i,bs) = gamma(i,bs,bs);
    end
end

% Calculates the BS that will transmit to the edge user
[~, bs_edge_user] = max(gamma_ib(N_users,:));

% % ========== Feasiblity test =============
% [P_ib, status] = min_power_global_cvx(Pt, BW, w, R, R_min_JT_user, gamma, isJT); % TODO: Do not take the max power into account?
% 
% % Initial power allocation
% P_ib = Pvec2mat(gamma, isJT, P_ib);
% P_ib_old = P_ib;

P_ib = zeros(N_users,N_BSs);
P_ib_old = zeros(N_users,N_BSs);
EE_total = [];
P_BS1 = [];
P_BS2 = [];
stop_BSs = false(N_BSs,1);
% if( status2exitflag(status) < 0)
%     l = 0;
% else
    for l = 0:L_max
        if(l~=0)
            p_aux = P_ib;
            p_aux(isnan(P_ib)) = 0;
            p_old_aux = P_ib_old;
            p_old_aux(isnan(P_ib_old)) = 0;
            %sum(abs(p_aux(:) - p_old_aux(:)).^2)
%             if(sum(abs(p_aux(:) - p_old_aux(:)).^2) < 1e-10 || EXITFLAG<=0)
%                 break;
%             end
            for bs = 1:N_BSs
                if(sum(abs(p_aux(:,bs) - p_old_aux(:,bs)).^2) < 1e-10 || EXITFLAG<=0)
                    stop_BSs(bs) = true;
                end
            end
            if(all(stop_BSs))
                break;
            end
            P_ib_old = P_ib;
        end
        
        for bs = 1:N_BSs
            if (~stop_BSs(bs))
                %SINR for non-CoMP case
                % TODO: add a parameter to say which BS transmits to the edge
                % user?
                %[ICI, INUI] = interference(gamma, isJT, P_ib(:));
                [ICI, ~] = interference_v2(gamma, isJT, P_ib);

                R_min_i = R(:,bs);
                if (bs == bs_edge_user)
                    R_min_i = [R_min_i; R_min_JT_user];
                end
                users_range = 1:(N_users - (bs~=bs_edge_user));
                [P_i, EXITFLAG, ~] = EE_NOMA_power_allocation_local_single_cell_cvx(Pt, BW, w, R_min_i, gamma_ib(users_range,bs), rho, P_fix, kappa, ICI(users_range,bs));
                if(EXITFLAG<0)
                    break;
                end
                %[EE, INUI_i, throughput_single_cell, p_tot_single_cell, user_i_data] = single_cell_EE(gamma_ib(users_range,bs), w, BW, P_i, ICI(users_range,bs), rho, kappa, P_fix);
                % Adds the P_i to the matrix of all BSs powers
                %sprintf("Total: %d | increment %d | BS %i iter %i",sum(P_i),sum(P_i) - nansum(P_ib(:,bs)),bs,l)

                %=== Calculate EE ===
                h_i = gamma_ib(users_range,bs);
                ICI_i = ICI(users_range,bs);

                % inter-NOMA-user interference (INUI)
                for k = 1:length(h_i)
                    INUI_i(k)=0;
                    for j = 1:(k-1)
                        INUI_i(k) = INUI_i(k) + P_i(j).*h_i(k);
                    end
                end 

                %Power consumption
                D = 0;
                for i = 1:length(h_i)
                    D = D + P_i(i) + P_i(i).*rho + kappa*(length(h_i)-i+1); 
                end
                D = D + P_fix;

                % Calculate SE
                N = 0;
                for i = 1:length(h_i)
                    N = N + w*BW*log(1+(h_i(i)*P_i(i))/(ICI_i(i) + INUI_i(i) + w));
                end

                EE = N/D;

                %=== End - Calculate EE ===


                if(length(P_i) ~= N_users)
                    P_ib(:,bs) = [P_i; NaN];
                else
                    P_ib(:,bs) = P_i;
                end
            end
        end
        
        %bs=1;
        %sprintf("TotalP: %d | BS %i iter %i", nansum(P_ib(:,bs)),bs,l)
        %bs=2;
        %sprintf("TotalP: %d | BS %i iter %i", nansum(P_ib(:,bs)),bs,l)
        % Calculate EE here
        P_i_aux = P_ib(:);
        P_i_aux(isnan(P_i_aux)) = 0;
        [R_tot,Ri] = system_throughput(w, BW, gamma, true, P_i_aux);
        PC = system_power_consumption(P_i_aux, gamma, rho, P_fix, kappa, true, false);
        
        EE_total_aux = R_tot/PC;
        x_axis = 1:l+1;
        EE_total = [EE_total; EE_total_aux];
        P_BS1 = [P_BS1; nansum(P_ib(:,1))];
        P_BS2 = [P_BS2; nansum(P_ib(:,2))];
        
%         %plot(x_axis,EE_total,'-+b')
%         %hold on,
%         plot(x_axis,P_BS1,'-*g')
%         hold on,
%         plot(x_axis,P_BS2,'-ok')
%         %legend('EE','Power BS_1','Power BS_2')
%         legend('Power BS_1','Power BS_2')
%         drawnow
    end
% end
P_i = P_ib(:);
P_i(isnan(P_i)) = 0;


end