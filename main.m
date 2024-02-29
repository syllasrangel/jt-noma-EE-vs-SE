clc;
clear all;
%close all;

%------- Simulation Parameters -------
N_inner_users = 2;
N_JT_users = 1;     % Keep it 1. (The simulation is not prepared for more than 1 JT user)
N_BSs = 2;
cell_radius = 0.6;
inner_circle_radius = 0.4;
h_length = 1;
P_dBm = 43;         % BS tranmit power in dBm (43 for macro cell)
N0_dBm = -139;      % Noise spectral density in dBm/Hz
BW_KHz = 180;       % Bandwidth per RB in KHz
%w = 5; % Total number of RBs (shared between BSs)
w = 100; % Total number of RBs (shared between BSs)
P_tol_dBm = 10;     % Minimum difference between the received power (with a normalized noise power) of a decoded signal and other non-decoded INUI signals
R_Kbps = [10, 500, 1000, 1500, 2000, 2500, 3000];
N_samples = 100;
d_u1 = 0.03;
PCM = "proposed"; % "M1", "M2", "proposed"
rho = 0.1;
kappa = 0.5; % 0, 0.5 and 2.5 
P_fix_dBm = 30; 

max_iter = 50; % Max number of iterations of the "ping-pong" algorithms
R_edge_div = 1; 

%Define the x axis parameter (some change might be needed in the code)
x_axis = R_Kbps;

%------- Conversions -------
Pt = (10^(P_dBm/10))/1000;          % [W]
N0 = (10^(N0_dBm/10))/1000;         % [W/Hz]
BW = BW_KHz*10^3;                   % [Hz]
P_tol = (10^(P_tol_dBm/10))/1000;   % [W]
R = R_Kbps*10^3;                    % [bps]
P_fix = (10^(P_fix_dBm/10))/1000;   % [W]

%------- Adjusts the PCM parameters -------
kappa_aux = kappa;
rho_aux = rho;
P_fix_aux = P_fix;
if PCM == "M1"
    rho = 0;
    P_fix = 0;
    kappa = 0;
end
if PCM == "M2"
    kappa = 0;
    rho = 0;
end

%------- Inicializations -------
R_EE_S1_global = zeros(length(x_axis),N_inner_users+N_JT_users,N_BSs,N_samples);
R_EE_S3_global = zeros(length(x_axis),N_inner_users+N_JT_users,N_BSs,N_samples);
R_EE_S1_global_SIC = zeros(length(x_axis),N_inner_users+N_JT_users,N_BSs,N_samples);
R_EE_S3_global_SIC = zeros(length(x_axis),N_inner_users+N_JT_users,N_BSs,N_samples);

Pi_EE_local = zeros(length(x_axis),2*(N_inner_users+N_JT_users),N_samples);
Exit_EE_local = zeros(length(x_axis),N_samples);
N_iter_local = zeros(length(x_axis),N_samples);
R_tot_EE_local = zeros(length(x_axis),N_samples);
R_EE_local = zeros(length(x_axis),N_inner_users+N_JT_users,N_BSs,N_samples);
PC_EE_local = zeros(length(x_axis),N_samples);

EE_total_iter =  zeros(length(x_axis),100,N_samples);
P_BS1_iter =  zeros(length(x_axis),100,N_samples);
P_BS2_iter =  zeros(length(x_axis),100,N_samples);

Pi_EE_S1_global_SIC = zeros(length(x_axis),2*(N_inner_users+N_JT_users),N_samples);
Pi_EE_S3_global_SIC = zeros(length(x_axis),2*N_inner_users+N_JT_users,N_samples);

Pi_EE_S1_global = zeros(length(x_axis),2*(N_inner_users+N_JT_users),N_samples);
Pi_EE_S3_global = zeros(length(x_axis),2*N_inner_users+N_JT_users,N_samples);

Exit_EE_S1_global = zeros(length(x_axis),N_samples);
Exit_EE_S3_global = zeros(length(x_axis),N_samples);

PC_EE_S1_global = zeros(length(x_axis),N_samples);
PC_EE_S3_global = zeros(length(x_axis),N_samples);
PC_EE_S1_global_SIC = zeros(length(x_axis),N_samples);
PC_EE_S3_global_SIC = zeros(length(x_axis),N_samples);

R_tot_EE_S1_global_1 = zeros(length(x_axis),N_samples);
R_tot_EE_S3_global_1 = zeros(length(x_axis),N_samples);
R_tot_EE_S1_global_1_SIC = zeros(length(x_axis),N_samples);
R_tot_EE_S3_global_1_SIC = zeros(length(x_axis),N_samples);

N_iter_S1_global = zeros(length(x_axis),N_samples);
N_iter_S3_global = zeros(length(x_axis),N_samples);

gamma_values = zeros(N_inner_users+N_JT_users,N_BSs,N_BSs,N_samples);

errors = {};

% Stores the date and time
date = datetime('now','TimeZone','local','Format','yyyy_MM_dd_HH_mm');

%for s=1:N_samples
s = 0;
N_errors = 0;

% Fixes the seed
rng(2);
error_samples = zeros(N_samples,1);

objective = "EE"; %"EE"; "SE"; "min"

while (s < N_samples)
    s = s + 1
    %------- Scenario -------
    % Scenario with 2 inner users in given positions
    %[BS1_position, BS2_position, inner_users_pos_BS1, inner_users_pos_BS2, JT_users] = generate_scenario_1(N_inner_users, N_JT_users, cell_radius, inner_circle_radius, [d_u1; 0.2]);
    
    % Scenario with N_inner_users inner users in random positions
    [BS1_position, BS2_position, inner_users_pos_BS1, inner_users_pos_BS2, JT_users] = generate_scenario_1(N_inner_users, N_JT_users, cell_radius, inner_circle_radius,[]);
    
    BSs_pos = [BS1_position; BS2_position];
    inner_users_pos = [inner_users_pos_BS1, inner_users_pos_BS2];
    JT_user_pos = JT_users;
    
    
    % Sort users according to their distance to their serving BS
    for bs = 1:N_BSs
        [inner_users_pos_dist,sort_idx] = sort(abs(BSs_pos(bs)-inner_users_pos(:,bs)), 'ascend');
        inner_users_pos(:,bs) = inner_users_pos(sort_idx,bs);
    end
    
    %------- Scenario Plot -------
%     plot_scenario(BSs_pos(1), inner_users_pos(:,1), JT_user_pos, BSs_pos(2), inner_users_pos(:,2), [], cell_radius, inner_circle_radius)
    
    %------- Wireless channel -------
    not_ordered = true; % Ensure that the generated channels are ordered.
    while(not_ordered)
        h = zeros(N_inner_users+1,N_BSs,N_BSs); % h(a,b,c) means the channel between BS c and the a-th user served by BS b
        for bs1 = 1:N_BSs
            for bs2 = 1:N_BSs
                for j = 1:N_inner_users
                    [h(j, bs1, bs2),~] = wireless_channel(h_length,abs(BSs_pos(bs2)-inner_users_pos(j,bs1)));
                end
            end
        end
        for bs = 1:N_BSs % CoMP user
            [h(N_inner_users+1,:, bs),~] = wireless_channel(h_length,abs(BSs_pos(bs)-JT_user_pos));
        end
        
        gamma = (abs(h).^2)./(N0*BW);
        
        for bs = 1:N_BSs
            for j = 1:N_inner_users
              not_ordered = not_ordered && gamma(j+1,bs,bs) < gamma(j,bs,bs);
            end
        end
        not_ordered = ~not_ordered;
    end
    
    gamma_values(:,:,:,s) = gamma;
   

    for x_axis_idx = 1:length(x_axis)
        %------- Power Allocation -------
        R_min = R(x_axis_idx)*ones(N_inner_users, N_BSs); % inner-users minimum data rate
        R_min_JT_user = R(x_axis_idx); % Edge-user minimum data rate
        
       try
            SIC_constraint=false; %If true, it adds a SIC constraint to the optimization
            if(objective == "EE")
            [Pi_EE_S1_global(x_axis_idx,:,s), Exit_EE_S1_global(x_axis_idx,s), N_iter_S1_global(x_axis_idx,s)] = EE_NOMA_power_allocation_global_cvx(Pt, BW, w, R_min, R_min_JT_user, gamma, rho, P_fix, kappa, PCM, true, SIC_constraint);
            end
            if(objective == "SE")
            [Pi_EE_S1_global(x_axis_idx,:,s), Exit_EE_S1_global(x_axis_idx,s), N_iter_S1_global(x_axis_idx,s)] = SE_NOMA_power_allocation_global_cvx(Pt, BW, w, R_min, R_min_JT_user, gamma, true, SIC_constraint);
            end
            if(objective == "min")
            [Pi_EE_S1_global(x_axis_idx,:,s), status] = min_power_global_cvx(Pt, BW, w, R_min, R_min_JT_user, gamma, true, SIC_constraint);
            Exit_EE_S1_global(x_axis_idx,s) = status2exitflag(status);
            end

            if(objective == "EE")
            [Pi_EE_S3_global(x_axis_idx,:,s), Exit_EE_S3_global(x_axis_idx,s), N_iter_S3_global(x_axis_idx,s)] = EE_NOMA_power_allocation_global_cvx(Pt, BW, w, R_min, R_min_JT_user, gamma, rho, P_fix, kappa, PCM, false, SIC_constraint);
            end
            if(objective == "SE")
            [Pi_EE_S3_global(x_axis_idx,:,s), Exit_EE_S3_global(x_axis_idx,s), N_iter_S3_global(x_axis_idx,s)] = SE_NOMA_power_allocation_global_cvx(Pt, BW, w, R_min, R_min_JT_user, gamma, false, SIC_constraint);
            end
            if(objective == "min")
            [Pi_EE_S3_global(x_axis_idx,:,s), status] = min_power_global_cvx(Pt, BW, w, R_min, R_min_JT_user, gamma, false, SIC_constraint);
            Exit_EE_S3_global(x_axis_idx,s) = status2exitflag(status);
            end
            catch ME
            % delete samples with error
            N_errors = N_errors+1;
            %             ME.stack.file
            %             ME.stack.line
            errors{N_errors} = sprintf('Sample: %d, Error id: %s', s, ME.identifier);            
            %s = s - 1;
            error_samples(s) = 1
            Exit_EE_S1_global(:,s) = NaN(size(x_axis));
            Exit_EE_S3_global(:,s) = NaN(size(x_axis));
            Exit_EE_local(:,s) = NaN(size(x_axis));
            break;
            end
        
        
        %------- Data Rate Calculation -------
        %===== Scenario 1 - JT-CoMP-NOMA =====
%         % Objective of maximizing the throughput
%         [R_BS1_S1(x_axis_idx,:,s), R_BS2_S1(x_axis_idx,:,s)] = S1_data_rate(w, BW, Pi_S1(x_axis_idx, 1:N_inner_users+N_JT_users,s), Pi_S1(x_axis_idx,N_inner_users+N_JT_users+1:2*(N_inner_users+N_JT_users),s), gamma_BS1, gamma_BS2, [gamma_BS2_1_inner; gamma_BS2_JT_CoMP], [gamma_BS1_2_inner; gamma_BS1_JT_CoMP], N_inner_users, N_JT_users);
%         [R_BS1_S1_global(x_axis_idx,:,s), R_BS2_S1_global(x_axis_idx,:,s)] = S1_data_rate(w, BW, Pi_S1_global(x_axis_idx, 1:N_inner_users+N_JT_users,s), Pi_S1_global(x_axis_idx,N_inner_users+N_JT_users+1:2*(N_inner_users+N_JT_users),s), gamma_BS1, gamma_BS2, [gamma_BS2_1_inner; gamma_BS2_JT_CoMP], [gamma_BS1_2_inner; gamma_BS1_JT_CoMP], N_inner_users, N_JT_users);
% 
        % Objective of maximizing the EE
        %[R_BS1_EE_S1(x_axis_idx,:,s), R_BS2_EE_S1(x_axis_idx,:,s)] = S1_data_rate(w, BW, Pi_EE_S1(x_axis_idx, 1:N_inner_users+N_JT_users,s), Pi_EE_S1(x_axis_idx,N_inner_users+N_JT_users+1:2*(N_inner_users+N_JT_users),s), gamma_BS1, gamma_BS2, [gamma_BS2_1_inner; gamma_BS2_JT_CoMP], [gamma_BS1_2_inner; gamma_BS1_JT_CoMP], N_inner_users, N_JT_users);
        is_cvx = false;
        [R_tot_EE_S1_global_1(x_axis_idx,s),R_EE_S1_global(x_axis_idx,:,:,s)] = system_throughput(w, BW, gamma, true, Pi_EE_S1_global(x_axis_idx,:,s).');
        PC_EE_S1_global(x_axis_idx,s) = system_power_consumption(Pi_EE_S1_global(x_axis_idx,:,s).', gamma, rho, P_fix, kappa, true, is_cvx);
%         %SE
%         [R_tot_SE_S1_global_1(x_axis_idx,s),R_SE_S1_global(x_axis_idx,:,:,s)] = system_throughput(w, BW, gamma, true, Pi_SE_S1_global(x_axis_idx,:,s).');
%         PC_SE_S1_global(x_axis_idx,s) = system_power_consumption(Pi_SE_S1_global(x_axis_idx,:,s).', gamma, rho, P_fix, kappa, true, is_cvx);
        

        
        %[R_BS1_EE_S1_global(x_axis_idx,:,s), R_BS2_EE_S1_global(x_axis_idx,:,s)] = S1_data_rate(w, BW, Pi_EE_S1_global(x_axis_idx, 1:N_inner_users+N_JT_users,s), Pi_EE_S1_global(x_axis_idx,N_inner_users+N_JT_users+1:2*(N_inner_users+N_JT_users),s), gamma);
%         
%         %===== Scenario 3 - Conventional NOMA =====
%         % Objective of maximizing the throughput
%         [R_BS1_S3(x_axis_idx,:,s), R_BS2_S3(x_axis_idx,:,s)] = S3_data_rate(w, BW, Pi_S3(x_axis_idx,1:N_inner_users+N_JT_users,s), Pi_S3(x_axis_idx,N_inner_users+N_JT_users+1:2*N_inner_users+N_JT_users,s), gamma_BS1, gamma_BS2, [gamma_BS2_1_inner; gamma_BS2_JT_CoMP], [gamma_BS1_2_inner; gamma_BS1_JT_CoMP], N_inner_users, N_JT_users);
%         [R_BS1_S3_global(x_axis_idx,:,s), R_BS2_S3_global(x_axis_idx,:,s)] = S3_data_rate(w, BW, Pi_S3_global(x_axis_idx,1:N_inner_users+N_JT_users,s), Pi_S3_global(x_axis_idx,N_inner_users+N_JT_users+1:2*N_inner_users+N_JT_users,s), gamma_BS1, gamma_BS2, [gamma_BS2_1_inner; gamma_BS2_JT_CoMP], [gamma_BS1_2_inner; gamma_BS1_JT_CoMP], N_inner_users, N_JT_users);
%
%         % Objective of maximizing the EE
%         [R_BS1_EE_S3(x_axis_idx,:,s), R_BS2_EE_S3(x_axis_idx,:,s)] = S3_data_rate(w, BW, Pi_EE_S3(x_axis_idx,1:N_inner_users+N_JT_users,s), Pi_EE_S3(x_axis_idx,N_inner_users+N_JT_users+1:2*N_inner_users+N_JT_users,s), gamma_BS1, gamma_BS2, [gamma_BS2_1_inner; gamma_BS2_JT_CoMP], [gamma_BS1_2_inner; gamma_BS1_JT_CoMP], N_inner_users, N_JT_users);
%         [R_BS1_EE_S3_global(x_axis_idx,:,s), R_BS2_EE_S3_global(x_axis_idx,:,s)] = S3_data_rate(w, BW, Pi_EE_S3_global(x_axis_idx,1:N_inner_users+N_JT_users,s), Pi_EE_S3_global(x_axis_idx,N_inner_users+N_JT_users+1:2*N_inner_users+N_JT_users,s), gamma_BS1, gamma_BS2, [gamma_BS2_1_inner; gamma_BS2_JT_CoMP], [gamma_BS1_2_inner; gamma_BS1_JT_CoMP], N_inner_users, N_JT_users);
        [R_tot_EE_S3_global_1(x_axis_idx,s),R_EE_S3_global(x_axis_idx,:,:,s)] = system_throughput(w, BW, gamma, false, Pi_EE_S3_global(x_axis_idx,:,s).');
        PC_EE_S3_global(x_axis_idx,s) = system_power_consumption(Pi_EE_S3_global(x_axis_idx,:,s).', gamma, rho, P_fix, kappa, false, is_cvx);
        
        [R_tot_EE_S3_global_1_SIC(x_axis_idx,s),R_EE_S3_global_SIC(x_axis_idx,:,:,s)] = system_throughput(w, BW, gamma, false, Pi_EE_S3_global_SIC(x_axis_idx,:,s).');
        PC_EE_S3_global_SIC(x_axis_idx,s) = system_power_consumption(Pi_EE_S3_global_SIC(x_axis_idx,:,s).', gamma, rho, P_fix, kappa, false, is_cvx);
        
        
    end
%     end
    filename = sprintf('workspaces/%s_workspace_%s_PCM_%s_kappa_%s_rho_%s_%dsamp.mat',date,objective,PCM, replace(sprintf("%0.2f", kappa),".","_"), replace(sprintf("%0.2f", rho),".","_"),N_samples);
    save(filename)
end

plot_outage = false;
plot_iter = false;
plot_throughput = true;
plot_EE = false;
plot_EC = false;
plot_individual_EE = false;

plot_results(filename, plot_outage, plot_iter, plot_throughput, plot_EE, plot_EC, plot_individual_EE)

%paper_figures()