function [] = paper_figures()


    data_EE = load('workspaces\2023_02_09_14_37_workspace_EE_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    data_SE = load('workspaces\2023_02_13_17_00_workspace_SE_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    %data_min = load('workspaces\2023_02_13_20_11_workspace_min_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    data_min = load('workspaces\2023_02_16_11_35_workspace_min_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    %data_min = load('workspaces/2023_02_14_11_03_workspace_min_PCM_proposed_kappa_0_50_rho_0_10_1000samp.mat');
    data_min_3_inner_users = load('workspaces/2023_02_14_20_46_workspace_min_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    data_min_4_inner_users = load('workspaces/2023_02_14_20_55_workspace_min_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    data_min_5_inner_users = load('workspaces/2023_02_14_21_02_workspace_min_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    
    % R = 1 Mbps
    data_EE_half_R = load('workspaces/2023_04_24_22_19_workspace_EE_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    data_EE_double_R = load('workspaces/2023_04_24_21_37_workspace_EE_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    data_EE_double_R_stop_per_BS = load('workspaces/2023_04_25_16_02_workspace_EE_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    data_EE_half_R_stop_per_BS = load('workspaces/2023_04_25_23_13_workspace_EE_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    x_axis = data_EE.x_axis;
    
    
    aux_kappa = 0.5;
    aux_rho = 0.1;
    aux_P_fix = 1;
    [R_JT_EE, R_NOMA_EE, ~, Pi_sys_JT_EE, Pi_sys_NOMA_EE, ~, feasible_JT_EE, feasible_NOMA_EE, ~] = prepare_data_box(data_EE, aux_kappa, aux_rho, aux_P_fix);
    [R_JT_SE, R_NOMA_SE, ~, Pi_sys_JT_SE, Pi_sys_NOMA_SE, ~, feasible_JT_SE, feasible_NOMA_SE, ~] = prepare_data_box(data_SE, aux_kappa, aux_rho, aux_P_fix);
    [R_JT_min, R_NOMA_min, ~, Pi_sys_JT_min, Pi_sys_NOMA_min, ~, feasible_JT_min, feasible_NOMA_min, ~] = prepare_data_box(data_min, aux_kappa, aux_rho, aux_P_fix);
    
    
    non_error_samples = ~data_EE.error_samples.' & ~data_SE.error_samples.' & ~data_min.error_samples.';
    feasible_all = feasible_JT_EE & feasible_NOMA_EE & feasible_JT_SE & feasible_NOMA_SE & feasible_JT_min & feasible_NOMA_min & non_error_samples;
    
    
    % -------------------------------
    % Iterations
    % -------------------------------
    [R_JT_EE_half_R, R_NOMA_EE_half_R, ~, Pi_sys_JT_EE_half_R, Pi_sys_NOMA_EE_half_R, ~, feasible_JT_EE_half_R, feasible_NOMA_EE_half_R, feasible_LOCAL_EE_half_R] = prepare_data_box(data_EE_half_R, aux_kappa, aux_rho, aux_P_fix);
    [R_JT_EE_double_R, R_NOMA_EE_double_R, ~, Pi_sys_JT_EE_double_R, Pi_sys_NOMA_EE_double_R, ~, feasible_JT_EE_double_R, feasible_NOMA_EE_double_R, feasible_LOCAL_EE_double_R] = prepare_data_box(data_EE_double_R, aux_kappa, aux_rho, aux_P_fix);
    [R_JT_EE_double_R_stop_per_BS, R_NOMA_EE_double_R_stop_per_BS, ~, Pi_sys_JT_EE_double_R_stop_per_BS, Pi_sys_NOMA_EE_double_R_stop_per_BS, ~, feasible_JT_EE_double_R_stop_per_BS, feasible_NOMA_EE_double_R_stop_per_BS, feasible_LOCAL_EE_double_R_stop_per_BS] = prepare_data_box(data_EE_double_R_stop_per_BS, aux_kappa, aux_rho, aux_P_fix);
    [R_JT_EE_half_R_stop_per_BS, R_NOMA_EE_half_R_stop_per_BS, ~, Pi_sys_JT_EE_half_R_stop_per_BS, Pi_sys_NOMA_EE_half_R_stop_per_BS, ~, feasible_JT_EE_half_R_stop_per_BS, feasible_NOMA_EE_half_R_stop_per_BS, feasible_LOCAL_EE_half_R_stop_per_BS] = prepare_data_box(data_EE_half_R_stop_per_BS, aux_kappa, aux_rho, aux_P_fix);
    
%     figure, plot(x_axis,mean(data_EE_half_R.N_iter_local(:,feasible_LOCAL_EE_half_R),2),'-+b','LineWidth',1),
%     hold on 
%     plot(x_axis,mean(data_EE_double_R.N_iter_local(:,feasible_LOCAL_EE_double_R),2),'-*k','LineWidth',1),
%     hold on 
%     plot(x_axis,mean(data_EE_double_R_stop_per_BS.N_iter_local(:,feasible_LOCAL_EE_double_R_stop_per_BS),2),'-*k','LineWidth',1),
%     legend('R/2', '2R','2R stop per BS')
%     xlabel('Minimum data rate requirement (Kbps)');
%     ylabel('Number of iterations');
    non_error_samples_local = ~data_EE_half_R.error_samples.' & ~data_EE_double_R.error_samples.' & ~data_EE_double_R_stop_per_BS.error_samples.' & ~data_EE_half_R_stop_per_BS.error_samples.';
    feasible_LOCAL_all = feasible_LOCAL_EE_half_R & feasible_LOCAL_EE_double_R & feasible_LOCAL_EE_double_R_stop_per_BS & feasible_LOCAL_EE_half_R_stop_per_BS & non_error_samples_local;
    N_iter_half = data_EE_half_R.N_iter_local(1,feasible_LOCAL_all);
    N_iter_double = data_EE_double_R.N_iter_local(1,feasible_LOCAL_all);
    N_iter_double_stop_per_BS = data_EE_double_R_stop_per_BS.N_iter_local(1,feasible_LOCAL_all);
    N_iter_half_stop_per_BS = data_EE_half_R_stop_per_BS.N_iter_local(1,feasible_LOCAL_all);
    
    figure
    boxplot([N_iter_half.',N_iter_double.', N_iter_double_stop_per_BS.', N_iter_half_stop_per_BS.'],'Labels',{'Half','Double', 'Double stop per BS', 'Half stop per BS'})
    %save('workspaces/figures_data.mat','Ri_EE','Ri_SE','Ri_min','-append');
    

    ss = 2;
    P_BS1_iter = data_EE_double_R.P_BS1_iter(1,1:50,ss);
    P_BS2_iter = data_EE_double_R.P_BS2_iter(1,1:50,ss);
    EE_total_iter = data_EE_double_R.EE_total_iter(1,1:50,ss);
    
    idxs = find(P_BS1_iter);
    figure, plot(1:length(idxs),P_BS1_iter(idxs),'-+b','LineWidth',1),
    hold on,
    plot(1:length(idxs),P_BS2_iter(idxs),'-ok','LineWidth',1),
    legend('BS 1', 'BS 2')
    xlabel('Number of iterations');
    ylabel('Total Power');

    figure, plot(1:length(idxs),EE_total_iter(idxs),'-+b','LineWidth',1),
%     hold on,
%     plot(1:length(idxs),ones(1,length(idxs)).*data_EE_double_R.R_tot_EE_S1_global_1(1,ss)./data_EE_double_R.PC_EE_S1_global(1,ss),'-ok','LineWidth',1),
%      hold on,
%     plot(1:length(idxs),ones(1,length(idxs)).*data_EE_double_R.R_tot_EE_local(1,ss)./data_EE_double_R.PC_EE_local(1,ss),'-*c','LineWidth',1),
    xlabel('Number of iterations');
    ylabel('EE');

    
    figure
    boxplot([N_iter_double_stop_per_BS.', N_iter_half_stop_per_BS.'],'Labels',{'Half','Double'})
    
    % -------------------------------
    % Is SIC satisfied? 
    % -------------------------------
%     data_sel = data_min;
%     ss = 0;
%     for s = find(feasible_all)
%         ss = ss + 1;
%         for x_axis_idx = 1:length(x_axis)
%             isJT = true;
%             if(isJT)
%                 Pi = data_sel.Pi_EE_S1_global(x_axis_idx,:,s);
%             else
%                 Pi = data_sel.Pi_EE_S3_global(x_axis_idx,:,s);
%             end
%             SIC_ok(x_axis_idx, ss) = is_SIC_satisfied(Pi, data_sel.BW, data_sel.w, data_sel.R_min, data_sel.R_min_JT_user, data_sel.gamma, isJT);
%             
%             if(~SIC_ok)
%                 it_went_wrong = 1;
%             end
%         end
%     end
    % -------------------------------
    % Outage all scenarios 
    % -------------------------------
    dt = data_min; % obj func does not matter (same outage)
    if(dt.s~=dt.N_samples)
        interval = 1:dt.s-1;
    else
        interval = 1:dt.s;
    end
    simu_samps = zeros(1,dt.N_samples);
    simu_samps(interval) = 1;
    
    % Ignores samples with error in the optimization process
    simu_samps = simu_samps & ~dt.error_samples.';
    
    outage_JT = dt.Exit_EE_S1_global(:,simu_samps & sum(isnan(dt.Exit_EE_S1_global),1)==0) < 0;
    outage_NOMA = dt.Exit_EE_S3_global(:,simu_samps & sum(isnan(dt.Exit_EE_S3_global),1)==0) < 0;
    save('workspaces/figures_data.mat','x_axis','outage_JT','outage_NOMA');
    
    figure, 
    plot(x_axis,mean(outage_NOMA,2),'-.ok','LineWidth',1),
    hold on,
    plot(x_axis,mean(outage_JT,2),'-+b','LineWidth',1),
    hold on,
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Outage ratio');
    xticks(x_axis);
    legend('NOMA', 'JTCN', 'Location', 'northwest')
    
    
    % -------------------------------
    % TP vs min data rate
    % -------------------------------
    rates = 1:length(x_axis);
    R_JT_EE_aux = R_JT_EE(rates,feasible_all);
    R_NOMA_EE_aux = R_NOMA_EE(rates,feasible_all);
    R_JT_SE_aux = R_JT_SE(rates,feasible_all);
    R_NOMA_SE_aux = R_NOMA_SE(rates,feasible_all);
    R_JT_min_aux = R_JT_min(rates,feasible_all);
    R_NOMA_min_aux = R_NOMA_min(rates,feasible_all);
    [R_EE_t_S1_global, R_EE_t_CI_S1_global] = mean_confidence_interval(R_JT_EE_aux);
    [R_EE_t_S3_global, R_EE_t_CI_S3_global] = mean_confidence_interval(R_NOMA_EE_aux);
    [R_SE_t_S1_global, R_SE_t_CI_S1_global] = mean_confidence_interval(R_JT_SE_aux);
    [R_SE_t_S3_global, R_SE_t_CI_S3_global] = mean_confidence_interval(R_NOMA_SE_aux);
    [R_min_t_S1_global, R_min_t_CI_S1_global] = mean_confidence_interval(R_JT_min_aux);
    [R_min_t_S3_global, R_min_t_CI_S3_global] = mean_confidence_interval(R_NOMA_min_aux);
    
    save('workspaces/figures_data.mat','R_JT_EE_aux','R_JT_SE_aux','R_JT_min_aux','R_NOMA_EE_aux','R_NOMA_SE_aux','R_NOMA_min_aux','-append');
    
    figure;
    hAx=axes;
    errorbar(x_axis,R_EE_t_S1_global, R_EE_t_CI_S1_global(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,R_EE_t_S3_global,R_EE_t_CI_S3_global(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,R_SE_t_S1_global, R_SE_t_CI_S1_global(:,2),'-+c','LineWidth',1)
    hold on,
    errorbar(x_axis,R_SE_t_S3_global,R_SE_t_CI_S3_global(:,2),'-.og','LineWidth',1),
    hold on,
    errorbar(x_axis,R_min_t_S1_global, R_min_t_CI_S1_global(:,2),'-+m','LineWidth',1)
    hold on,
    errorbar(x_axis,R_min_t_S3_global,R_min_t_CI_S3_global(:,2),'-.oy','LineWidth',1),
    hold on,
    legend('EE JTCN','EE NOMA','SE JTCN','SE NOMA', 'min JTCN','min NOMA');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Network throughput (b/s)');
    %title(sprintf('Throughput - %s',PCM_text));
    hAx.YScale='log';
    
    % -------------------------------
    % EE vs min data rate
    % -------------------------------
    Pi_sys_JT_EE_aux = Pi_sys_JT_EE(rates,feasible_all);
    Pi_sys_NOMA_EE_aux = Pi_sys_NOMA_EE(rates,feasible_all);
    Pi_sys_JT_SE_aux = Pi_sys_JT_SE(rates,feasible_all);
    Pi_sys_NOMA_SE_aux = Pi_sys_NOMA_SE(rates,feasible_all);
    Pi_sys_JT_min_aux = Pi_sys_JT_min(rates,feasible_all);
    Pi_sys_NOMA_min_aux = Pi_sys_NOMA_min(rates,feasible_all);
    [EE_EE_t_S1_global, EE_EE_t_CI_S1_global] = mean_confidence_interval(R_JT_EE_aux./Pi_sys_JT_EE_aux);
    [EE_EE_t_S3_global, EE_EE_t_CI_S3_global] = mean_confidence_interval(R_NOMA_EE_aux./Pi_sys_NOMA_EE_aux);
    [EE_SE_t_S1_global, EE_SE_t_CI_S1_global] = mean_confidence_interval(R_JT_SE_aux./Pi_sys_JT_SE_aux);
    [EE_SE_t_S3_global, EE_SE_t_CI_S3_global] = mean_confidence_interval(R_NOMA_SE_aux./Pi_sys_NOMA_SE_aux);
    [EE_min_t_S1_global, EE_min_t_CI_S1_global] = mean_confidence_interval(R_JT_min_aux./Pi_sys_JT_min_aux);
    [EE_min_t_S3_global, EE_min_t_CI_S3_global] = mean_confidence_interval(R_NOMA_min_aux./Pi_sys_NOMA_min_aux);
    
    save('workspaces/figures_data.mat','Pi_sys_JT_EE_aux','Pi_sys_JT_SE_aux','Pi_sys_JT_min_aux','Pi_sys_NOMA_EE_aux','Pi_sys_NOMA_SE_aux','Pi_sys_NOMA_min_aux','-append');
    
    figure;
    hAx=axes;
    errorbar(x_axis,EE_EE_t_S1_global, EE_EE_t_CI_S1_global(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,EE_EE_t_S3_global,EE_EE_t_CI_S3_global(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,EE_SE_t_S1_global, EE_SE_t_CI_S1_global(:,2),'-+c','LineWidth',1)
    hold on,
    errorbar(x_axis,EE_SE_t_S3_global,EE_SE_t_CI_S3_global(:,2),'-.og','LineWidth',1),
    hold on,
    errorbar(x_axis,EE_min_t_S1_global, EE_min_t_CI_S1_global(:,2),'-+m','LineWidth',1)
    hold on,
    errorbar(x_axis,EE_min_t_S3_global,EE_min_t_CI_S3_global(:,2),'-.oy','LineWidth',1),
    hold on,
    legend('EE JTCN','EE NOMA','SE JTCN','SE NOMA', 'min JTCN','min NOMA');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Network energy efficiency (b/s/Joule)');
    %title(sprintf('Throughput - %s',PCM_text));
    hAx.YScale='log';
    
    % -------------------------------
    % System power consumption
    % -------------------------------
    [EC_EE_t_S1_global, EC_EE_t_CI_S1_global] = mean_confidence_interval(Pi_sys_JT_EE_aux);
    [EC_EE_t_S3_global, EC_EE_t_CI_S3_global] = mean_confidence_interval(Pi_sys_NOMA_EE_aux);
    [EC_SE_t_S1_global, EC_SE_t_CI_S1_global] = mean_confidence_interval(Pi_sys_JT_SE_aux);
    [EC_SE_t_S3_global, EC_SE_t_CI_S3_global] = mean_confidence_interval(Pi_sys_NOMA_SE_aux);
    [EC_min_t_S1_global, EC_min_t_CI_S1_global] = mean_confidence_interval(Pi_sys_JT_min_aux);
    [EC_min_t_S3_global, EC_min_t_CI_S3_global] = mean_confidence_interval(Pi_sys_NOMA_min_aux);
    
    figure;
    hAx=axes;
    errorbar(x_axis,EC_EE_t_S1_global, EC_EE_t_CI_S1_global(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,EC_EE_t_S3_global,EC_EE_t_CI_S3_global(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,EC_SE_t_S1_global, EC_SE_t_CI_S1_global(:,2),'-+c','LineWidth',1)
    hold on,
    errorbar(x_axis,EC_SE_t_S3_global,EC_SE_t_CI_S3_global(:,2),'-.og','LineWidth',1),
    hold on,
    errorbar(x_axis,EC_min_t_S1_global, EC_min_t_CI_S1_global(:,2),'-+m','LineWidth',1)
    hold on,
    errorbar(x_axis,EC_min_t_S3_global,EC_min_t_CI_S3_global(:,2),'-.oy','LineWidth',1),
    hold on,
    legend('EE JTCN','EE NOMA','SE JTCN','SE NOMA', 'min JTCN','min NOMA');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Average system power consumption');
    %title(sprintf('Throughput - %s',PCM_text));
    %hAx.YScale='log';
    
    % -------------------------------
    % Boxplot TP
    % -------------------------------
    rate = 4;
    R_JT_EE_aux = R_JT_EE(rate,feasible_all);
    R_NOMA_EE_aux = R_NOMA_EE(rate,feasible_all);
    R_JT_SE_aux = R_JT_SE(rate,feasible_all);
    R_NOMA_SE_aux = R_NOMA_SE(rate,feasible_all);
    R_JT_min_aux = R_JT_min(rate,feasible_all);
    R_NOMA_min_aux = R_NOMA_min(rate,feasible_all);
    %save('workspaces/figures_data.mat','R_JT_050','R_NOMA_050','R_LOCAL_050','-append');
    
    %box_data = {(R_JT_050(rates,:)).', (R_NOMA_050(rates,:)).', (R_LOCAL_050(rates,:)).'};
    
    box_data = {[R_JT_EE_aux; R_JT_SE_aux; R_JT_min_aux].', [R_NOMA_EE_aux; R_NOMA_SE_aux; R_NOMA_min_aux].'};
    figure,
    boxplotGroup(box_data, 'PrimaryLabels', {'JTCN' 'NOMA'}, ...
      'SecondaryLabels',cellstr(["EE", "SE", "min"]), 'InterGroupSpace', 1, ...
       'GroupType','withinGroups')
    %title('Throughput');
    ylabel('Network throughput (b/s)')
    
    % -------------------------------
    % Boxplot EE
    % -------------------------------
    rate = 4;
    EE_JT_EE_aux = R_JT_EE(rate,feasible_all)./Pi_sys_JT_EE(rate,feasible_all);
    EE_NOMA_EE_aux = R_NOMA_EE(rate,feasible_all)./Pi_sys_NOMA_EE(rate,feasible_all);
    EE_JT_SE_aux = R_JT_SE(rate,feasible_all)./Pi_sys_JT_SE(rate,feasible_all);
    EE_NOMA_SE_aux = R_NOMA_SE(rate,feasible_all)./Pi_sys_NOMA_SE(rate,feasible_all);
    EE_JT_min_aux = R_JT_min(rate,feasible_all)./Pi_sys_JT_min(rate,feasible_all);
    EE_NOMA_min_aux = R_NOMA_min(rate,feasible_all)./Pi_sys_NOMA_min(rate,feasible_all);
    %save('workspaces/figures_data.mat','R_JT_050','R_NOMA_050','R_LOCAL_050','-append');
    
    box_data = {[EE_JT_EE_aux; EE_JT_SE_aux; EE_JT_min_aux].', [EE_NOMA_EE_aux; EE_NOMA_SE_aux; EE_NOMA_min_aux].'};
    figure,
    boxplotGroup(box_data, 'PrimaryLabels', {'JTCN' 'NOMA'}, ...
      'SecondaryLabels',cellstr(["EE", "SE", "min"]), 'InterGroupSpace', 1, ...
       'GroupType','withinGroups')
    %title('Throughput');
    ylabel('Network energy efficiency (b/s/Joule)')
    
    
    % -------------------------------
    % EE vs SE - gain in EE and loss in TP
    % -------------------------------    
    EE_EE = R_JT_EE(:,feasible_all)./Pi_sys_JT_EE(:,feasible_all);
    EE_SE = R_JT_SE(:,feasible_all)./Pi_sys_JT_SE(:,feasible_all);
    EE_gain = EE_EE./EE_SE;
    %TP_loss = R_JT_SE(:,feasible_all)./R_JT_EE(:,feasible_all);
    TP_gain = R_JT_EE(:,feasible_all)./R_JT_SE(:,feasible_all);
    %EC_decrease = Pi_sys_JT_SE(:,feasible_all)./Pi_sys_JT_EE(:,feasible_all);
    EC_increase = Pi_sys_JT_EE(:,feasible_all)./Pi_sys_JT_SE(:,feasible_all);
    
    [EE_gain_t, EE_gain_CI] = mean_confidence_interval(EE_gain);
    [TP_gain_t, TP_gain_CI] = mean_confidence_interval(TP_gain);
    [EC_increase_t, EC_increase_CI] = mean_confidence_interval(EC_increase);
    
    lin2dB = @(x)10*log10(x);
    
    figure;
    hAx=axes;
    errorbar(x_axis,lin2dB(EE_gain_t), EE_gain_CI(:,2),'-+g','LineWidth',1)
    hold on,
    errorbar(x_axis,lin2dB(TP_gain_t), TP_gain_CI(:,2),'-.or','LineWidth',1),
    hold on,
    errorbar(x_axis,lin2dB(EC_increase_t), EC_increase_CI(:,2),'-.ob','LineWidth',1),
    hold on,
    xticks(x_axis);
    legend('EE gain','TP gain','Power increase');
    xlabel('Minimum data rate requirement (Kbps)');
    %hAx.YScale='log';
    
    EE_gain_EE_over_SE = lin2dB(EE_gain_t);
    TP_gain_EE_over_SE = lin2dB(TP_gain_t);
    EC_gain_EE_over_SE = lin2dB(EC_increase_t);
    
    save('workspaces/figures_data.mat','EE_gain_EE_over_SE','TP_gain_EE_over_SE','EC_gain_EE_over_SE','-append');
    
    % -------------------------------
    % EE vs minP - gain in EE and loss in TP
    % -------------------------------  

    EE_EE = R_JT_EE(:,feasible_all)./Pi_sys_JT_EE(:,feasible_all);
    %EE_SE = R_JT_SE(:,feasible_all)./Pi_sys_JT_SE(:,feasible_all);
    EE_min = R_JT_min(:,feasible_all)./Pi_sys_JT_min(:,feasible_all);
    EE_gain = EE_EE./EE_min;
    TP_gain = R_JT_EE(:,feasible_all)./R_JT_min(:,feasible_all);
    EC_increase = Pi_sys_JT_EE(:,feasible_all)./Pi_sys_JT_min(:,feasible_all);
    %TP_loss = R_JT_SE(:,feasible_all)./R_JT_EE(:,feasible_all);
    %EC_decrease = Pi_sys_JT_SE(:,feasible_all)./Pi_sys_JT_EE(:,feasible_all);
    
    [EE_gain_t, EE_gain_CI] = mean_confidence_interval(EE_gain);
    [TP_gain_t, TP_gain_CI] = mean_confidence_interval(TP_gain);
    [EC_increase_t, EC_increase_CI] = mean_confidence_interval(EC_increase);
    
    figure;
    hAx=axes;
%     errorbar(x_axis,lin2dB(EE_gain_t), EE_gain_CI(:,2),'-+g','LineWidth',1)
%     hold on,
%     errorbar(x_axis,lin2dB(TP_gain_t), TP_gain_CI(:,2),'-.or','LineWidth',1),
%     hold on,
%     errorbar(x_axis,lin2dB(EC_increase_t), EC_increase_CI(:,2),'-.ob','LineWidth',1),
%     hold on,
    plot(x_axis,lin2dB(EE_gain_t),'-+g','LineWidth',1)
    hold on,
    plot(x_axis,lin2dB(TP_gain_t),'-.or','LineWidth',1),
    hold on,
    plot(x_axis,lin2dB(EC_increase_t),'-.ob','LineWidth',1),
    hold on,
    xticks(x_axis);
    legend('EE gain','TP gain','Power increase');
    xlabel('Minimum data rate requirement (Kbps)');
    %hAx.YScale='log';
    
    EE_gain_EE_over_min = lin2dB(EE_gain_t);
    TP_gain_EE_over_min = lin2dB(TP_gain_t);
    EC_gain_EE_over_min = lin2dB(EC_increase_t);
    
    save('workspaces/figures_data.mat','EE_gain_EE_over_min','TP_gain_EE_over_min','EC_gain_EE_over_min','-append');
    
    % -------------------------------
    % SumRate vs minP - gain in EE and loss in TP
    % -------------------------------  

    %EE_EE = R_JT_EE(:,feasible_all)./Pi_sys_JT_EE(:,feasible_all);
    EE_SE = R_JT_SE(:,feasible_all)./Pi_sys_JT_SE(:,feasible_all);
    EE_min = R_JT_min(:,feasible_all)./Pi_sys_JT_min(:,feasible_all);
    EE_gain = EE_SE./EE_min;
    TP_gain = R_JT_SE(:,feasible_all)./R_JT_min(:,feasible_all);
    EC_increase = Pi_sys_JT_SE(:,feasible_all)./Pi_sys_JT_min(:,feasible_all);
    %TP_loss = R_JT_SE(:,feasible_all)./R_JT_EE(:,feasible_all);
    %EC_decrease = Pi_sys_JT_SE(:,feasible_all)./Pi_sys_JT_EE(:,feasible_all);
    
    [EE_gain_t, EE_gain_CI] = mean_confidence_interval(EE_gain);
    [TP_gain_t, TP_gain_CI] = mean_confidence_interval(TP_gain);
    [EC_increase_t, EC_increase_CI] = mean_confidence_interval(EC_increase);
    
    figure;
    hAx=axes;
    plot(x_axis,lin2dB(EE_gain_t),'-+g','LineWidth',1)
    hold on,
    plot(x_axis,lin2dB(TP_gain_t),'-.or','LineWidth',1),
    hold on,
    plot(x_axis,lin2dB(EC_increase_t),'-.ob','LineWidth',1),
    hold on,
    xticks(x_axis);
    legend('EE gain','TP gain','Power increase');
    xlabel('Minimum data rate requirement (Kbps)');
    %hAx.YScale='log';
    
    EE_gain_SE_over_min = lin2dB(EE_gain_t);
    TP_gain_SE_over_min = lin2dB(TP_gain_t);
    EC_gain_SE_over_min = lin2dB(EC_increase_t);
    
    save('workspaces/figures_data.mat','EE_gain_SE_over_min','TP_gain_SE_over_min','EC_gain_SE_over_min','-append');
    
    % -------------------------------
    % Fairness
    % -------------------------------   
    R_ib_EE=data_EE.R_EE_S1_global;
    R_ib_SE=data_SE.R_EE_S1_global;
    R_ib_min=data_min.R_EE_S1_global;
    R_ib_EE(isnan(R_ib_EE)) = 0;
    R_ib_SE(isnan(R_ib_SE)) = 0;
    R_ib_min(isnan(R_ib_min)) = 0;
    N_users = data_EE.N_inner_users*data_EE.N_BSs + 1;
    jains_index_EE = (squeeze(sum(sum(R_ib_EE,2),3)).^2)./(N_users * squeeze(sum(sum(R_ib_EE.^2,2),3)));
    jains_index_SE = (squeeze(sum(sum(R_ib_SE,2),3)).^2)./(N_users * squeeze(sum(sum(R_ib_SE.^2,2),3)));
    jains_index_min = (squeeze(sum(sum(R_ib_min,2),3)).^2)./(N_users * squeeze(sum(sum(R_ib_min.^2,2),3)));
    
    jains_index_EE=jains_index_EE(:,feasible_all);
    jains_index_SE=jains_index_SE(:,feasible_all);
    jains_index_min=jains_index_min(:,feasible_all);
    save('workspaces/figures_data.mat','jains_index_EE','jains_index_SE','jains_index_min','-append');
    
    [jains_index_EE_t, jains_index_EE_t_CI] = mean_confidence_interval(jains_index_EE);
    [jains_index_SE_t, jains_index_SE_t_CI] = mean_confidence_interval(jains_index_SE);
    [jains_index_min_t, jains_index_min_t_CI] = mean_confidence_interval(jains_index_min);    
    
    figure;
    hAx=axes;
    errorbar(x_axis,jains_index_EE_t, jains_index_EE_t_CI(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,jains_index_SE_t,jains_index_SE_t_CI(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,jains_index_min_t, jains_index_min_t_CI(:,2),'-+c','LineWidth',1)
    hold on,
    legend('EE','TP','minP');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel("Average jain's fairness index");
    
    % -------------------------------
    % Users data rate - Boxplot
    % -------------------------------
    R_idx = 4;
    Ri_EE = [squeeze(R_ib_EE(R_idx,:,1, feasible_all)); squeeze(R_ib_EE(R_idx,1:data_EE.N_inner_users,2, feasible_all))];
    Ri_SE = [squeeze(R_ib_SE(R_idx,:,1, feasible_all)); squeeze(R_ib_SE(R_idx,1:data_SE.N_inner_users,2, feasible_all))];
    Ri_min = [squeeze(R_ib_min(R_idx,:,1, feasible_all)); squeeze(R_ib_min(R_idx,1:data_min.N_inner_users,2, feasible_all))];

    save('workspaces/figures_data.mat','Ri_EE','Ri_SE','Ri_min','-append');
    
    box_data = {Ri_EE.', Ri_SE.', Ri_min.'};
    figure,boxplotGroup(box_data, 'PrimaryLabels', {'EE' 'TP', 'min'}, ...
      'SecondaryLabels',cellstr(["User 1","User 2", "edge user", "User 1 (BS2)","User 2 (BS2)"]), 'InterGroupSpace', 1, ...
       'GroupType','withinGroups')
    ylabel('User data rate (b/s)')
    
    
    
 