function [] = paper_figures()


    data_PCM1 = load("workspaces\2022_12_08_18_31_workspace_PCM_M1_kappa_0_00_rho_0_00_100samp.mat");
    data_PCM2 = load("workspaces\2022_12_09_00_25_workspace_PCM_M2_kappa_0_00_rho_0_00_100samp.mat");
    %data_kappa_00 = load("workspaces\2022_12_12_15_50_workspace_PCM_proposed_kappa_0_00_rho_0_10_100samp.mat"); % Without local
    data_kappa_00 = load("workspaces\2022_12_13_21_37_workspace_PCM_proposed_kappa_0_00_rho_0_10_100samp.mat");
    data_kappa_050 = load('workspaces\2022_12_08_18_32_workspace_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    data_kappa_250 = load("workspaces\2022_12_11_23_18_workspace_PCM_proposed_kappa_2_50_rho_0_10_100samp.mat");
    
    x_axis = data_PCM1.x_axis;


    aux_kappa = 0.5;
    aux_rho = 0.1;
    aux_P_fix = 1;
    [R_JT_PCM1, R_NOMA_PCM1, R_LOCAL_PCM1, Pi_sys_JT_PCM1, Pi_sys_NOMA_PCM1, Pi_sys_LOCALl_PCM1, feasible_JT_PCM1, feasible_NOMA_PCM1, feasible_LOCAL_PCM1] = prepare_data_box(data_PCM1, aux_kappa, aux_rho, aux_P_fix);
    [R_JT_PCM2, R_NOMA_PCM2, R_LOCAL_PCM2, Pi_sys_JT_PCM2, Pi_sys_NOMA_PCM2, Pi_sys_LOCAL_PCM2, feasible_JT_PCM2, feasible_NOMA_PCM2, feasible_LOCAL_PCM2] = prepare_data_box(data_PCM2, aux_kappa, aux_rho, aux_P_fix);
    [R_JT_050, R_NOMA_050, R_LOCAL_050, Pi_sys_JT_050, Pi_sys_NOMA_050, Pi_sys_LOCAL_050, feasible_JT_050, feasible_NOMA_050, feasible_LOCAL_050] = prepare_data_box(data_kappa_050, aux_kappa, aux_rho, aux_P_fix);
    %[R_JT_250, R_NOMA_250, R_LOCAL_250, Pi_sys_JT_250, Pi_sys_NOMA_250, Pi_sys_local_250, feasible_JT_250, feasible_NOMA_250, feasible_LOCAL_250] = prepare_data_box(data_kappa_250, aux_kappa, aux_rho, aux_P_fix);
    
    feasible = feasible_JT_PCM1 & feasible_JT_PCM2 & feasible_JT_050; %& feasible_JT_250;
      
    % -------------------------------
    % EE opt. PCM 1 eval all PCMs
    % -------------------------------
    [~, ~, ~, Pi_sys_JT_eval_PCM1, ~, ~, ~, ~, ~] = prepare_data_box(data_PCM1, 0, 0, 0);
    [~, ~, ~, Pi_sys_JT_eval_PCM2, ~, ~, ~, ~, ~] = prepare_data_box(data_PCM1, 0, 0, aux_P_fix);
    [~, ~, ~, Pi_sys_JT_eval_PCMprop, ~, ~, ~, ~, ~] = prepare_data_box(data_PCM1, aux_kappa, aux_rho, aux_P_fix);
    [EE_EE_t_S1_global_PCM_prop, EE_EE_t_CI_S1_global_PCM_prop] = mean_confidence_interval(R_JT_PCM1(:,feasible_JT_PCM1)./Pi_sys_JT_eval_PCMprop(:,feasible_JT_PCM1));
    [EE_EE_t_S1_global_PCM2, EE_EE_t_CI_S1_global_PCM2] = mean_confidence_interval(R_JT_PCM1(:,feasible_JT_PCM1)./Pi_sys_JT_eval_PCM2(:,feasible_JT_PCM1));
    [EE_EE_t_S1_global_PCM1, EE_EE_t_CI_S1_global_PCM1] = mean_confidence_interval(R_JT_PCM1(:,feasible_JT_PCM1)./Pi_sys_JT_eval_PCM1(:,feasible_JT_PCM1));
    
    y_EE_opt_1_eval_050 = R_JT_PCM1(:,feasible_JT_PCM1)./Pi_sys_JT_eval_PCMprop(:,feasible_JT_PCM1);
    y_EE_opt_1_eval_2 = R_JT_PCM1(:,feasible_JT_PCM1)./Pi_sys_JT_eval_PCM2(:,feasible_JT_PCM1);
    y_EE_opt_1_eval_1 = R_JT_PCM1(:,feasible_JT_PCM1)./Pi_sys_JT_eval_PCM1(:,feasible_JT_PCM1);
    
    save('workspaces/figures_data.mat','x_axis','y_EE_opt_1_eval_050','y_EE_opt_1_eval_2','y_EE_opt_1_eval_1');
    
    figure;
    hAx=axes;
    errorbar(x_axis,EE_EE_t_S1_global_PCM1, EE_EE_t_CI_S1_global_PCM1(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,EE_EE_t_S1_global_PCM2, EE_EE_t_CI_S1_global_PCM2(:,2),'-.ok','LineWidth',1)
    hold on,
    errorbar(x_axis,EE_EE_t_S1_global_PCM_prop, EE_EE_t_CI_S1_global_PCM_prop(:,2),'-.og','LineWidth',1)
    hold on,
    legend('Eval. PCM 1', 'Eval. PCM 2', 'Eval. PCM Proposed');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Average energy efficiency (b/s/Joule)');
    xticks(x_axis);
    title(sprintf('EE JT-CoMP-NOMA - Optimized with PCM-1 (rho = %0.2f, kappa = %0.2f, P_{fix} = %0.2f)', aux_rho, aux_kappa, aux_P_fix));
    hAx.YScale='log';
    
    sprintf('EE PCM-1/PCMprop for %d Kbps = %0.2f',x_axis(1), EE_EE_t_S1_global_PCM1(1)/EE_EE_t_S1_global_PCM_prop(1))
    sprintf('EE PCM-1/PCMprop for %d Kbps = %0.2f',x_axis(5), EE_EE_t_S1_global_PCM1(5)/EE_EE_t_S1_global_PCM_prop(5))
    sprintf('EE PCM-1/PCM-2 for %d Kbps = %0.2f',x_axis(1), EE_EE_t_S1_global_PCM1(1)/EE_EE_t_S1_global_PCM2(1))
    sprintf('EE PCM-1/PCM-2 for %d Kbps = %0.2f',x_axis(5), EE_EE_t_S1_global_PCM1(5)/EE_EE_t_S1_global_PCM2(5))   
    % -------------------------------
    % TP opt. all
    % -------------------------------
    [TP_EE_t_JT_PCM1, TP_EE_t_CI_JT_PCM1] = mean_confidence_interval(R_JT_PCM1(:,feasible));
    [TP_EE_t_JT_PCM2, TP_EE_t_CI_JT_PCM2] = mean_confidence_interval(R_JT_PCM2(:,feasible));
    [TP_EE_t_JT_prop, TP_EE_t_CI_JT_prop] = mean_confidence_interval(R_JT_050(:,feasible));
    %[EE_EE_t_JT_prop_250, EE_EE_t_CI_JT_prop_250] = mean_confidence_interval(R_JT_250(:,feasible));
    
    y_TP_opt_050_eval_050 = R_JT_050(:,feasible);
    y_TP_opt_2_eval_050 = R_JT_PCM2(:,feasible);
    y_TP_opt_1_eval_050 = R_JT_PCM1(:,feasible);
    
    
    minTP = (data_PCM1.N_BSs*data_PCM1.N_inner_users + data_PCM1.N_JT_users)*data_PCM1.R;
    
    save('workspaces/figures_data.mat','y_TP_opt_050_eval_050','y_TP_opt_2_eval_050','y_TP_opt_1_eval_050','minTP','-append');
    
    
    figure,
    hAx=axes;
    errorbar(x_axis,TP_EE_t_JT_PCM1, TP_EE_t_CI_JT_PCM1(:,2),'-+b','LineWidth',1),
    hold on,
    errorbar(x_axis,TP_EE_t_JT_PCM2, TP_EE_t_CI_JT_PCM2(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,TP_EE_t_JT_prop, TP_EE_t_CI_JT_prop(:,2),'-.og','LineWidth',1),
    hold on,
    plot(x_axis, minTP,'-r','LineWidth',1)
    legend('Opt. PCM 1','Opt. PCM 2','Opt. PCM prop. \kappa = 0.50','Location', 'southeast');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Average throughput (b/s)');
    xticks(x_axis);
    title(sprintf('JT - Evaluated with PCM prop (rho = %0.2f, kappa = %0.2f, P_{fix} = %0.2f)', aux_rho, aux_kappa, aux_P_fix));
    hAx.YScale='log';
    
    aloc_power_JT_PCM1 = squeeze(sum(data_PCM1.Pi_EE_S1_global,2));
    aloc_power_JT_PCM2 = squeeze(sum(data_PCM2.Pi_EE_S1_global,2));    
    aloc_power_JT_PCMprop = squeeze(sum(data_kappa_050.Pi_EE_S1_global,2));
            
    [avg_aloc_power_JT_PCM1, ~] = mean_confidence_interval(aloc_power_JT_PCM1(:, feasible_JT_PCM1));
    [avg_aloc_power_JT_PCM2, ~] = mean_confidence_interval(aloc_power_JT_PCM2(:, feasible_JT_PCM2));
    [avg_aloc_power_JT_PCMprop, ~] = mean_confidence_interval(aloc_power_JT_PCMprop(:, feasible_JT_050));
    
    sprintf('TP PCMprop/PCM-1 for %d Kbps = %0.2f',x_axis(1), TP_EE_t_JT_prop(1)/TP_EE_t_JT_PCM1(1))
    sprintf('TP PCMprop/PCM-1 for %d Kbps = %0.2f',x_axis(7), TP_EE_t_JT_prop(7)/TP_EE_t_JT_PCM1(7))
    
    sprintf('Avg. aloc power PCM-1 for %d Kbps = %0.2f (%0.2f%%)',x_axis(1), avg_aloc_power_JT_PCM1(1), 100*avg_aloc_power_JT_PCM1(1)/(data_PCM1.N_BSs*data_PCM1.Pt))
    sprintf('Avg. aloc power PCMprop for %d Kbps = %0.2f (%0.2f%%)',x_axis(1), avg_aloc_power_JT_PCMprop(1), 100*avg_aloc_power_JT_PCMprop(1)/(data_PCM1.N_BSs*data_PCM1.Pt))
    sprintf('Avg. aloc power PCM-1 for %d Kbps = %0.2f (%0.2f%%)',x_axis(7), avg_aloc_power_JT_PCM1(7), 100*avg_aloc_power_JT_PCM1(7)/(data_PCM1.N_BSs*data_PCM1.Pt))
    sprintf('Avg. aloc power PCMprop for %d Kbps = %0.2f (%0.2f%%)',x_axis(7), avg_aloc_power_JT_PCMprop(7), 100*avg_aloc_power_JT_PCMprop(7)/(data_PCM1.N_BSs*data_PCM1.Pt))
    
     % -------------------------------
    % EE opt. all eval PCM prop
    % -------------------------------
    [EE_EE_t_JT_PCM1, EE_EE_t_CI_JT_PCM1] = mean_confidence_interval(R_JT_PCM1(:,feasible)./Pi_sys_JT_PCM1(:,feasible));
    [EE_EE_t_JT_PCM2, EE_EE_t_CI_JT_PCM2] = mean_confidence_interval(R_JT_PCM2(:,feasible)./Pi_sys_JT_PCM2(:,feasible));
    [EE_EE_t_JT_prop, EE_EE_t_CI_JT_prop] = mean_confidence_interval(R_JT_050(:,feasible)./Pi_sys_JT_050(:,feasible));
    %[EE_EE_t_JT_prop_250, EE_EE_t_CI_JT_prop_250] = mean_confidence_interval(R_JT_250(:,feasible)./Pi_sys_JT_250(:,feasible));
    
    y_EE_opt_050_eval_050 = R_JT_050(:,feasible)./Pi_sys_JT_050(:,feasible);
    y_EE_opt_2_eval_050 = R_JT_PCM2(:,feasible)./Pi_sys_JT_PCM2(:,feasible);
    y_EE_opt_1_eval_050 = R_JT_PCM1(:,feasible)./Pi_sys_JT_PCM1(:,feasible);
    
    save('workspaces/figures_data.mat','y_EE_opt_050_eval_050','y_EE_opt_2_eval_050','y_EE_opt_1_eval_050','-append');
    
    figure,
    hAx=axes;
    errorbar(x_axis,EE_EE_t_JT_PCM1, EE_EE_t_CI_JT_PCM1(:,2),'-+b','LineWidth',1),
    hold on,
    errorbar(x_axis,EE_EE_t_JT_PCM2, EE_EE_t_CI_JT_PCM2(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,EE_EE_t_JT_prop, EE_EE_t_CI_JT_prop(:,2),'-.og','LineWidth',1),
    hold on,
    legend('Opt. PCM 1','Opt. PCM 2','Opt. PCM prop. \kappa = 0.50');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Average energy efficiency (b/s/Joule)');
    xticks(x_axis);
    title(sprintf('JT - Evaluated with PCM prop (rho = %0.2f, kappa = %0.2f, P_{fix} = %0.2f)', aux_rho, aux_kappa, aux_P_fix));
    %hAx.YScale='log';
    
    sprintf('EE PCMprop/PCM-1 for %d Kbps = %0.2f',x_axis(1), EE_EE_t_JT_prop(1)/EE_EE_t_JT_PCM1(1))
    sprintf('EE PCMprop/PCM-1 for %d Kbps = %0.2f',x_axis(7), EE_EE_t_JT_prop(7)/EE_EE_t_JT_PCM1(7))
    sprintf('EE PCMprop/PCM-2 for %d Kbps = %0.2f',x_axis(1), EE_EE_t_JT_prop(1)/EE_EE_t_JT_PCM2(1))
    sprintf('EE PCMprop/PCM-2 for %d Kbps = %0.2f',x_axis(7), EE_EE_t_JT_prop(7)/EE_EE_t_JT_PCM2(7))
    
    % -------------------------------
    % Outage all scenarios 
    % -------------------------------
    dt = data_PCM2; % PCM does not matter (same outage)
    if(dt.s~=dt.N_samples)
        interval = 1:dt.s-1;
    else
        interval = 1:dt.s;
    end
    simu_samps = zeros(1,dt.N_samples);
    simu_samps(interval) = 1;
    
    outage_JT = dt.Exit_EE_S1_global(:,simu_samps & sum(isnan(dt.Exit_EE_S1_global),1)==0) < 0;
    outage_NOMA = dt.Exit_EE_S3_global(:,simu_samps & sum(isnan(dt.Exit_EE_S3_global),1)==0) < 0;
    outage_LOCAL = dt.Exit_EE_local(:,simu_samps & sum(isnan(dt.Exit_EE_local),1)==0) < 0;
    save('workspaces/figures_data.mat','outage_JT','outage_NOMA','outage_LOCAL','-append');
    
    figure, 
    plot(x_axis,mean(outage_NOMA,2),'-.ok','LineWidth',1),
    hold on,
    plot(x_axis,mean(outage_JT,2),'-+b','LineWidth',1),
    hold on,
    plot(x_axis,mean(outage_LOCAL,2),'--og','LineWidth',1)    
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Outage ratio');
    xticks(x_axis);
    legend('Conventional NOMA', 'JT-CoMP NOMA', 'Local DPS-CoMP NOMA','Location', 'northwest')
    
    
    % -------------------------------
    % Boxplot EE opt. and eval with PCM prop
    % -------------------------------
    rates = [1,4,7];
    feasible_all = feasible_JT_050 & feasible_NOMA_050 & feasible_LOCAL_050;
    EE_JT_050 = (R_JT_050(:,feasible_all)./Pi_sys_JT_050(:,feasible_all));
    EE_NOMA_050 = (R_NOMA_050(:,feasible_all)./Pi_sys_NOMA_050(:,feasible_all));
    EE_LOCAL_050 = (R_LOCAL_050(:,feasible_all)./Pi_sys_LOCAL_050(:,feasible_all));
    
    save('workspaces/figures_data.mat','EE_JT_050','EE_NOMA_050','EE_LOCAL_050','-append');
    
    box_data = {EE_JT_050(rates,:).', EE_NOMA_050(rates,:).', EE_LOCAL_050(rates,:).'};
    figure,boxplotGroup(box_data, 'PrimaryLabels', {'JTCN' 'Conv.', 'ILO'}, ...
      'SecondaryLabels',cellstr(string(x_axis(rates))), 'InterGroupSpace', 1, ...
       'GroupType','withinGroups')
    title('Energy Efficiency');
    ylabel('Average energy efficiency (b/s/Joule)')
    
    
    sprintf('EE JTCN/NOMA for %d Kbps = %0.5f',x_axis(1), mean(R_JT_050(1,feasible_all)./Pi_sys_JT_050(1,feasible_all))./mean(R_NOMA_050(1,feasible_all)./Pi_sys_NOMA_050(1,feasible_all)))
    sprintf('EE JTCN/NOMA for %d Kbps = %0.2f',x_axis(4), mean(R_JT_050(4,feasible_all)./Pi_sys_JT_050(4,feasible_all))./mean(R_NOMA_050(4,feasible_all)./Pi_sys_NOMA_050(4,feasible_all)))
    sprintf('EE JTCN/NOMA for %d Kbps = %0.2f',x_axis(7), mean(R_JT_050(7,feasible_all)./Pi_sys_JT_050(7,feasible_all))./mean(R_NOMA_050(7,feasible_all)./Pi_sys_NOMA_050(7,feasible_all)))
    sprintf('EE JTCN/ILO for %d Kbps = %0.5f',x_axis(1), mean(R_JT_050(1,feasible_all)./Pi_sys_JT_050(1,feasible_all))./mean(R_LOCAL_050(1,feasible_all)./Pi_sys_LOCAL_050(1,feasible_all)))
    sprintf('EE JTCN/ILO for %d Kbps = %0.2f',x_axis(4), mean(R_JT_050(4,feasible_all)./Pi_sys_JT_050(4,feasible_all))./mean(R_LOCAL_050(4,feasible_all)./Pi_sys_LOCAL_050(4,feasible_all)))
    sprintf('EE JTCN/ILO for %d Kbps = %0.2f',x_axis(7), mean(R_JT_050(7,feasible_all)./Pi_sys_JT_050(7,feasible_all))./mean(R_LOCAL_050(7,feasible_all)./Pi_sys_LOCAL_050(7,feasible_all)))
    % -------------------------------
    % Boxplot TP opt. with PCM prop
    % -------------------------------
    rates = [1,4,7];
    R_JT_050 = R_JT_050(:,feasible_all);
    R_NOMA_050 = R_NOMA_050(:,feasible_all);
    R_LOCAL_050 = R_LOCAL_050(:,feasible_all);
    save('workspaces/figures_data.mat','R_JT_050','R_NOMA_050','R_LOCAL_050','-append');
    
    box_data = {(R_JT_050(rates,:)).', (R_NOMA_050(rates,:)).', (R_LOCAL_050(rates,:)).'};
    figure,boxplotGroup(box_data, 'PrimaryLabels', {'JT' 'Conv.', 'ILO'}, ...
      'SecondaryLabels',cellstr(string(x_axis(rates))), 'InterGroupSpace', 1, ...
       'GroupType','withinGroups')
    %title('Throughput');
    ylabel('Average throughput (b/s)')

    
    % -------------------------------
    % Boxplot EE var. kappa
    % -------------------------------
    [R_JT_00, R_NOMA_00, R_LOCAL_00, Pi_sys_JT_00, Pi_sys_NOMA_00, Pi_sys_local_00, feasible_JT_00, feasible_NOMA_00, feasible_LOCAL_00] = prepare_data_box(data_kappa_00);
    [R_JT_050, R_NOMA_050, R_LOCAL_050, Pi_sys_JT_050, Pi_sys_NOMA_050, Pi_sys_local_050, feasible_JT_050, feasible_NOMA_050, feasible_LOCAL_050] = prepare_data_box(data_kappa_050);
    [R_JT_250, R_NOMA_250, R_LOCAL_250, Pi_sys_JT_250, Pi_sys_NOMA_250, Pi_sys_local_250, feasible_JT_250, feasible_NOMA_250, feasible_LOCAL_250] = prepare_data_box(data_kappa_250);
    
    R_idx = 4;
    feasible = feasible_JT_00 & feasible_JT_050 & feasible_JT_250 & feasible_NOMA_00 & feasible_NOMA_050 & feasible_NOMA_250;
    JT_00 = (R_JT_00(R_idx,feasible)./Pi_sys_JT_00(R_idx,feasible)).';
    JT_050 = (R_JT_050(R_idx,feasible)./Pi_sys_JT_050(R_idx,feasible)).';
    JT_250 = (R_JT_250(R_idx,feasible)./Pi_sys_JT_250(R_idx,feasible)).';
    NOMA_00 = (R_NOMA_00(R_idx,feasible)./Pi_sys_NOMA_00(R_idx,feasible)).';
    NOMA_050 = (R_NOMA_050(R_idx,feasible)./Pi_sys_NOMA_050(R_idx,feasible)).';
    NOMA_250 = (R_NOMA_250(R_idx,feasible)./Pi_sys_NOMA_250(R_idx,feasible)).';
    
    EE_JT_diff_kappa = [JT_00, JT_050, JT_250].';
    EE_NOMA_diff_kappa = [NOMA_00, NOMA_050, NOMA_250].';
    save('workspaces/figures_data.mat','EE_JT_diff_kappa','EE_NOMA_diff_kappa','-append');
    
    box_data = {[JT_00, JT_050, JT_250], [NOMA_00, NOMA_050, NOMA_250]};
    figure,boxplotGroup(box_data, 'PrimaryLabels', {'JT' 'Conv.'}, ...
      'SecondaryLabels',cellstr(["\kappa = 0","\kappa = 0.5","\kappa = 2.5"]), 'InterGroupSpace', 1, ...
       'GroupType','withinGroups')
    %title('Energy Efficiency');
    ylabel('Average energy efficiency (b/s/Joule)')

    sprintf('EE JTCN/NOMA for %d Kbps and kappa 0 = %0.2f',x_axis(R_idx), mean(JT_00)./mean(NOMA_00))
    sprintf('EE JTCN/NOMA for %d Kbps and kappa 0.5 = %0.2f',x_axis(R_idx), mean(JT_050)./mean(NOMA_050))
    sprintf('EE JTCN/NOMA for %d Kbps and kappa 2.5 = %0.2f',x_axis(R_idx), mean(JT_250)./mean(NOMA_250))
    
    aloc_power_JT_00 = squeeze(sum(data_kappa_00.Pi_EE_S1_global,2));
    aloc_power_JT_050 = squeeze(sum(data_kappa_050.Pi_EE_S1_global,2));    
    aloc_power_JT_250 = squeeze(sum(data_kappa_250.Pi_EE_S1_global,2));
     
    fea = feasible_JT_00 && feasible_NOMA_00;
    [avg_aloc_power_JT_00, ~] = mean_confidence_interval(aloc_power_JT_00(:, feasible_NOMA_00));
    [avg_aloc_power_JT_050, ~] = mean_confidence_interval(aloc_power_JT_050(:, feasible_NOMA_050));
    [avg_aloc_power_JT_250, ~] = mean_confidence_interval(aloc_power_JT_250(:, feasible_NOMA_250));
    
    [avg_Pi_sys_JT_00, ~] = mean_confidence_interval(Pi_sys_JT_00(:, feasible_NOMA_00));
    [avg_Pi_sys_JT_050, ~] = mean_confidence_interval(Pi_sys_JT_050(:, feasible_NOMA_050));
    [avg_Pi_sys_JT_250, ~] = mean_confidence_interval(Pi_sys_JT_250(:, feasible_NOMA_250));
    
    user_PC_JT_00 = 0.1*avg_aloc_power_JT_00;
    SIC_PC_00 =  0;
    PC_circuit = 2*1;
    
    avg_aloc_power_JT_00 + user_PC_JT_00 + SIC_PC_00 + PC_circuit
    
    avg_Pi_sys_JT_00
    
    user_PC_JT_050 = 0.1*avg_aloc_power_JT_050;
    SIC_PC_050 =  0.5*9;
    PC_circuit = 2*1;
    
    avg_aloc_power_JT_050 + user_PC_JT_050 + SIC_PC_050 + PC_circuit
    
    user_PC_JT_250 = 0.1*avg_aloc_power_JT_250;
    SIC_PC_250 =  2.5*9;
    PC_circuit = 2*1;
    
    avg_aloc_power_JT_250 + user_PC_JT_250 + SIC_PC_250 + PC_circuit
    
    
    
    % -------------------------------
    % Users data rate - Boxplot
    % -------------------------------
    R_EE_S1_global = data_kappa_050.R_EE_S1_global;
    R_EE_S3_global = data_kappa_050.R_EE_S3_global;
    R_EE_S1_global(isnan(R_EE_S1_global)) = 0;
    R_EE_S3_global(isnan(R_EE_S3_global)) = 0;

    R_idx = 4;
    feasible = feasible_JT_050 & feasible_NOMA_050;
    Ri_JT = R_EE_S1_global(R_idx,:,1, feasible);
    Ri_NOMA = R_EE_S3_global(R_idx,:,1, feasible);
    %Ri_JT = [squeeze(Ri_JT).', squeeze(R_EE_S1_global(R_idx,1:data_kappa_050.N_inner_users,2, feasible)).'];
    %Ri_NOMA = [squeeze(Ri_NOMA).', squeeze(R_EE_S3_global(R_idx,1:data_kappa_050.N_inner_users,2, feasible)).'];
    Ri_JT = [squeeze(Ri_JT); squeeze(R_EE_S1_global(R_idx,1:data_kappa_050.N_inner_users,2, feasible))];
    Ri_NOMA = [squeeze(Ri_NOMA); squeeze(R_EE_S3_global(R_idx,1:data_kappa_050.N_inner_users,2, feasible))];
    
    save('workspaces/figures_data.mat','Ri_JT','Ri_NOMA','-append');
    
    box_data = {Ri_JT.', Ri_NOMA.'};
    figure,boxplotGroup(box_data, 'PrimaryLabels', {'JT' 'Conv.'}, ...
      'SecondaryLabels',cellstr(["User 1","User 2", "edge user", "User 1 (BS2)","User 2 (BS2)"]), 'InterGroupSpace', 1, ...
       'GroupType','withinGroups')
    ylabel('User data rate (b/s)')

    % -------------------------------
    % Users EE
    % -------------------------------
    %TODO: set the same yaxis limit for all figs
    data = data_kappa_050;
    N_users = data.N_inner_users+data.N_JT_users;
    
    feasible = feasible_JT_00 & feasible_NOMA_00 & feasible_JT_050 & feasible_NOMA_050 & feasible_JT_250 & feasible_NOMA_250;
    
    tot_user_power_i_bs_JT_00 = zeros(length(x_axis),N_users,data_kappa_00.N_BSs,data_kappa_00.N_samples);
    tot_user_power_i_bs_NOMA_00 = zeros(length(x_axis),N_users,data_kappa_00.N_BSs,data_kappa_00.N_samples);
    tot_user_power_i_bs_JT_050 = zeros(length(x_axis),N_users,data_kappa_050.N_BSs,data_kappa_050.N_samples);
    tot_user_power_i_bs_NOMA_050 = zeros(length(x_axis),N_users,data_kappa_050.N_BSs,data_kappa_050.N_samples);
    tot_user_power_i_bs_JT_250 = zeros(length(x_axis),N_users,data_kappa_250.N_BSs,data_kappa_250.N_samples);
    tot_user_power_i_bs_NOMA_250 = zeros(length(x_axis),N_users,data_kappa_250.N_BSs,data_kappa_250.N_samples);
%     tot_user_power_i_bs_S1 = zeros(length(x_axis),N_users,data.N_BSs,data.N_samples);
%     tot_user_power_i_bs_S3 = zeros(length(x_axis),N_users,data.N_BSs,data.N_samples);
    for samp = 1:data.N_samples
        for x_idx = 1:length(x_axis)
            for bs=1:data.N_BSs
%                 Pib_EE_S1_global = Pvec2mat(data.gamma, true, data.Pi_EE_S1_global(x_idx,:,samp));
%                 Pib_EE_S3_global = Pvec2mat(data.gamma, false, data.Pi_EE_S3_global(x_idx,:,samp));
%                 [~, tot_user_power_i_bs_S1(x_idx,:,bs,samp)] = user_power_consumption_2(Pib_EE_S1_global(:,bs), data.gamma, data.rho, data.kappa, false);
%                 [~, tot_user_power_i_bs_S3(x_idx,:,bs,samp)] = user_power_consumption_2(Pib_EE_S3_global(:,bs), data.gamma, data.rho, data.kappa, false);
                
                Pib_JT_00 = Pvec2mat(data_kappa_00.gamma, true, data_kappa_00.Pi_EE_S1_global(x_idx,:,samp));
                Pib_NOMA_00 = Pvec2mat(data_kappa_00.gamma, false, data_kappa_00.Pi_EE_S3_global(x_idx,:,samp));
                [~, tot_user_power_i_bs_JT_00(x_idx,:,bs,samp)] = user_power_consumption_2(Pib_JT_00(:,bs), data_kappa_00.gamma, data_kappa_00.rho, data_kappa_00.kappa, false);
                [~, tot_user_power_i_bs_NOMA_00(x_idx,:,bs,samp)] = user_power_consumption_2(Pib_NOMA_00(:,bs), data_kappa_00.gamma, data_kappa_00.rho, data_kappa_00.kappa, false);
                
                Pib_JT_050 = Pvec2mat(data_kappa_050.gamma, true, data_kappa_050.Pi_EE_S1_global(x_idx,:,samp));
                Pib_NOMA_050 = Pvec2mat(data_kappa_050.gamma, false, data_kappa_050.Pi_EE_S3_global(x_idx,:,samp));
                [~, tot_user_power_i_bs_JT_050(x_idx,:,bs,samp)] = user_power_consumption_2(Pib_JT_050(:,bs), data_kappa_050.gamma, data_kappa_050.rho, data_kappa_050.kappa, false);
                [~, tot_user_power_i_bs_NOMA_050(x_idx,:,bs,samp)] = user_power_consumption_2(Pib_NOMA_050(:,bs), data_kappa_050.gamma, data_kappa_050.rho, data_kappa_050.kappa, false);
                
                Pib_JT_250 = Pvec2mat(data_kappa_250.gamma, true, data_kappa_250.Pi_EE_S1_global(x_idx,:,samp));
                Pib_NOMA_250 = Pvec2mat(data_kappa_250.gamma, false, data_kappa_250.Pi_EE_S3_global(x_idx,:,samp));
                [~, tot_user_power_i_bs_JT_250(x_idx,:,bs,samp)] = user_power_consumption_2(Pib_JT_250(:,bs), data_kappa_250.gamma, data_kappa_250.rho, data_kappa_250.kappa, false);
                [~, tot_user_power_i_bs_NOMA_250(x_idx,:,bs,samp)] = user_power_consumption_2(Pib_NOMA_250(:,bs), data_kappa_250.gamma, data_kappa_250.rho, data_kappa_250.kappa, false);
            end
        end
    end
    
    R_JT_00 = data_kappa_00.R_EE_S1_global;
    R_NOMA_00 = data_kappa_00.R_EE_S3_global;
    R_JT_050 = data_kappa_050.R_EE_S1_global;
    R_NOMA_050 = data_kappa_050.R_EE_S3_global;
    R_JT_250 = data_kappa_250.R_EE_S1_global;
    R_NOMA_250 = data_kappa_250.R_EE_S3_global;
    R_JT_00(isnan(R_JT_00)) = 0;
    R_NOMA_00(isnan(R_NOMA_00)) = 0;
    R_JT_050(isnan(R_JT_050)) = 0;
    R_NOMA_050(isnan(R_NOMA_050)) = 0;
    R_JT_250(isnan(R_JT_250)) = 0;
    R_NOMA_250(isnan(R_NOMA_250)) = 0;
    % For each cell-center user
    EEi_JT_00 = zeros(data.N_BSs*data.N_inner_users + 1, length(x_axis), sum(feasible));
    EEi_NOMA_00 = zeros(data.N_BSs*data.N_inner_users + 1, length(x_axis), sum(feasible));
    EEi_JT_050 = zeros(data.N_BSs*data.N_inner_users + 1, length(x_axis), sum(feasible));
    EEi_NOMA_050 = zeros(data.N_BSs*data.N_inner_users + 1, length(x_axis), sum(feasible));
    EEi_JT_250 = zeros(data.N_BSs*data.N_inner_users + 1, length(x_axis), sum(feasible));
    EEi_NOMA_250 = zeros(data.N_BSs*data.N_inner_users + 1, length(x_axis), sum(feasible));
    for bs=1:data.N_BSs
        for ii = 1:N_users-1
            EEi_JT_00(two_dim_2_one_dim(ii, bs, N_users, false), :, :) = squeeze(R_JT_00(:,ii,bs,feasible)./tot_user_power_i_bs_JT_00(:,ii,bs,feasible));
            EEi_NOMA_00(two_dim_2_one_dim(ii, bs, N_users, false), :, :) = squeeze(R_NOMA_00(:,ii,bs,feasible)./tot_user_power_i_bs_NOMA_00(:,ii,bs,feasible));
            EEi_JT_050(two_dim_2_one_dim(ii, bs, N_users, false), :, :) = squeeze(R_JT_050(:,ii,bs,feasible)./tot_user_power_i_bs_JT_050(:,ii,bs,feasible));
            EEi_NOMA_050(two_dim_2_one_dim(ii, bs, N_users, false), :, :) = squeeze(R_NOMA_050(:,ii,bs,feasible)./tot_user_power_i_bs_NOMA_050(:,ii,bs,feasible));
            EEi_JT_250(two_dim_2_one_dim(ii, bs, N_users, false), :, :) = squeeze(R_JT_250(:,ii,bs,feasible)./tot_user_power_i_bs_JT_250(:,ii,bs,feasible));
            EEi_NOMA_250(two_dim_2_one_dim(ii, bs, N_users, false), :, :) = squeeze(R_NOMA_250(:,ii,bs,feasible)./tot_user_power_i_bs_NOMA_250(:,ii,bs,feasible));
            
            [EE_JT_00, EE_JT_00_CI] = mean_confidence_interval( squeeze(R_JT_00(:,ii,bs,feasible)./tot_user_power_i_bs_JT_00(:,ii,bs,feasible)) );
            [EE_NOMA_00, EE_NOMA_00_CI] = mean_confidence_interval( squeeze(R_NOMA_00(:,ii,bs,feasible)./tot_user_power_i_bs_NOMA_00(:,ii,bs,feasible)) );
            [EE_JT_050, EE_JT_050_CI] = mean_confidence_interval( squeeze(R_JT_050(:,ii,bs,feasible)./tot_user_power_i_bs_JT_050(:,ii,bs,feasible)) );
            [EE_NOMA_050, EE_NOMA_050_CI] = mean_confidence_interval( squeeze(R_NOMA_050(:,ii,bs,feasible)./tot_user_power_i_bs_NOMA_050(:,ii,bs,feasible)) );
            [EE_JT_250, EE_JT_250_CI] = mean_confidence_interval( squeeze(R_JT_250(:,ii,bs,feasible)./tot_user_power_i_bs_JT_250(:,ii,bs,feasible)) );
            [EE_NOMA_250, EE_NOMA_250_CI] = mean_confidence_interval( squeeze(R_NOMA_250(:,ii,bs,feasible)./tot_user_power_i_bs_NOMA_250(:,ii,bs,feasible)) );
            
            
            figure;
            hAx=axes;
            errorbar(x_axis,EE_JT_00,EE_JT_00_CI(:,2),'-+b','LineWidth',1)
            hold on,
            errorbar(x_axis,EE_NOMA_00,EE_NOMA_00_CI(:,2),'-.ok','LineWidth',1),
            hold on,
            errorbar(x_axis,EE_JT_050,EE_JT_050_CI(:,2),'-+b','LineWidth',1)
            hold on,
            errorbar(x_axis,EE_NOMA_050,EE_NOMA_050_CI(:,2),'-.ok','LineWidth',1),
            hold on,
            errorbar(x_axis,EE_JT_250,EE_JT_250_CI(:,2),'-+b','LineWidth',1)
            hold on,
            errorbar(x_axis,EE_NOMA_250,EE_NOMA_250_CI(:,2),'-.ok','LineWidth',1),
            hold on,
            legend('EE JT-CoMP-NOMA','EE Conventional NOMA');
            xlabel('Minimum data rate requirement (Kbps)');
            ylabel('Average energy efficiency (b/s/Joule)');
            hAx.YScale='log';
            title(sprintf('User %d BS %d',ii,bs));
        end
    end
    
    [EE_JT_00, EE_JT_00_CI] = mean_confidence_interval( squeeze(R_JT_00(:,N_users,1,feasible)./sum(tot_user_power_i_bs_JT_00(:,N_users,:,feasible),3)) );
    [EE_NOMA_00, EE_NOMA_00_CI] = mean_confidence_interval( squeeze(R_NOMA_00(:,N_users,1,feasible)./tot_user_power_i_bs_NOMA_00(:,N_users,1,feasible)) );
    [EE_JT_050, EE_JT_050_CI] = mean_confidence_interval( squeeze(R_JT_050(:,N_users,1,feasible)./sum(tot_user_power_i_bs_JT_050(:,N_users,:,feasible),3)) );
    [EE_NOMA_050, EE_NOMA_050_CI] = mean_confidence_interval( squeeze(R_NOMA_050(:,N_users,1,feasible)./tot_user_power_i_bs_NOMA_050(:,N_users,1,feasible)) );
    [EE_JT_250, EE_JT_250_CI] = mean_confidence_interval( squeeze(R_JT_250(:,N_users,1,feasible)./sum(tot_user_power_i_bs_JT_250(:,N_users,:,feasible),3)) );
    [EE_NOMA_250, EE_NOMA_250_CI] = mean_confidence_interval( squeeze(R_NOMA_250(:,N_users,1,feasible)./tot_user_power_i_bs_NOMA_250(:,N_users,1,feasible)) );
    
    EEi_JT_00(two_dim_2_one_dim(N_users, 1, N_users, false), :, :) = squeeze(R_JT_00(:,N_users,1,feasible)./sum(tot_user_power_i_bs_JT_00(:,N_users,:,feasible),3));
    EEi_NOMA_00(two_dim_2_one_dim(N_users, 1, N_users, false), :, :) = squeeze(R_NOMA_00(:,N_users,1,feasible)./tot_user_power_i_bs_NOMA_00(:,N_users,1,feasible));
    EEi_JT_050(two_dim_2_one_dim(N_users, 1, N_users, false), :, :) = squeeze(R_JT_050(:,N_users,1,feasible)./sum(tot_user_power_i_bs_JT_050(:,N_users,:,feasible),3));
    EEi_NOMA_050(two_dim_2_one_dim(N_users, 1, N_users, false), :, :) = squeeze(R_NOMA_050(:,N_users,1,feasible)./tot_user_power_i_bs_NOMA_050(:,N_users,1,feasible));
    EEi_JT_250(two_dim_2_one_dim(N_users, 1, N_users, false), :, :) = squeeze(R_JT_250(:,N_users,1,feasible)./sum(tot_user_power_i_bs_JT_250(:,N_users,:,feasible),3));
    EEi_NOMA_250(two_dim_2_one_dim(N_users, 1, N_users, false), :, :) = squeeze(R_NOMA_250(:,N_users,1,feasible)./tot_user_power_i_bs_NOMA_250(:,N_users,1,feasible));
    
    save('workspaces/figures_data.mat','EEi_JT_00','EEi_NOMA_00','EEi_JT_050','EEi_NOMA_050','EEi_JT_250','EEi_NOMA_250','-append');
    
    figure;
    hAx=axes;
    errorbar(x_axis,EE_JT_00,EE_JT_00_CI(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,EE_NOMA_00,EE_NOMA_00_CI(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,EE_JT_050,EE_JT_050_CI(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,EE_NOMA_050,EE_NOMA_050_CI(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,EE_JT_250,EE_JT_250_CI(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,EE_NOMA_250,EE_NOMA_250_CI(:,2),'-.ok','LineWidth',1),
    hold on,
    legend('EE JT-CoMP-NOMA','EE Conventional NOMA');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Average energy efficiency (b/s/Joule)');
    hAx.YScale='log';
    title(sprintf('User %d (Edge User)',N_users)); 
    
end