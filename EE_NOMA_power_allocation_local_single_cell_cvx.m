function [P_i, EXITFLAG, l] = EE_NOMA_power_allocation_local_single_cell_cvx(Pt, BW, w, R_min_i, h_i, rho, P_fix, kappa, ICI_i)
L_max=30; % Max number of iterations

N_users = length(h_i);

% Sets the inital power as the minimum required to attend the QoS requirement
% TODO: check calculations and when to use R_min_JT_user
P_i_old = zeros(N_users,1);
for i = 1:N_users
    INUI = 0;
    for j=1:i-1
        INUI = INUI + P_i_old(j)*h_i(i);
    end
    P_i_old(i) = (2^(R_min_i(i)/(w*BW))-1).*(ICI_i(i) + INUI + w)./h_i(i);
end
P_i = P_i_old;

for l = 0:L_max

    if(l~=0)
        if(abs(single_cell_EE(h_i, w, BW, P_i_old, ICI_i, rho, kappa, P_fix) - single_cell_EE(h_i, w, BW, P_i, ICI_i, rho, kappa, P_fix)) < 0.0001 || status2exitflag(status)<0)
            break;
        end
        P_i_old = P_i;
    end

    [~, INUI_i] = single_cell_EE(h_i, w, BW, P_i, ICI_i, rho, kappa, P_fix);

    SINR_i = P_i.*h_i./(ICI_i + INUI_i + w);
    a_i = SINR_i./(1+SINR_i);
    c_i = log2(1+SINR_i) - SINR_i.*log2(SINR_i)./(1+SINR_i);
    c_i(isnan(c_i))=0;

    [P_i, status, ~] = dinkelbach_algorithm_local_cvx(Pt, BW, w, R_min_i, h_i, rho, P_fix, kappa, a_i, c_i, ICI_i);
end

EXITFLAG = status2exitflag(status);
end