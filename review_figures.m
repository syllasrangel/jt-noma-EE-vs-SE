function [] = review_figures()

    % R = 1 Mbps
%     data_EE_half_R = load('workspaces/2023_04_24_22_19_workspace_EE_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
%     data_EE_double_R = load('workspaces/2023_04_24_21_37_workspace_EE_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
%     data_EE_double_R_stop_per_BS = load('workspaces/2023_04_25_16_02_workspace_EE_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
%     data_EE_half_R_stop_per_BS = load('workspaces/2023_04_25_23_13_workspace_EE_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');

    data_EE_half_R = load('workspaces/2023_04_25_23_13_workspace_EE_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    data_EE_double_R = load('workspaces/2023_04_25_16_02_workspace_EE_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    data_EE = load('workspaces/2023_04_25_23_57_workspace_EE_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
    
    aux_kappa = 0.5;
    aux_rho = 0.1;
    aux_P_fix = 1;


    % -------------------------------
    % Iterations
    % -------------------------------
    [R_JT_EE_half_R, R_NOMA_EE_half_R, ~, Pi_sys_JT_EE_half_R, Pi_sys_NOMA_EE_half_R, ~, feasible_JT_EE_half_R, feasible_NOMA_EE_half_R, feasible_LOCAL_EE_half_R] = prepare_data_box(data_EE_half_R, aux_kappa, aux_rho, aux_P_fix);
    [R_JT_EE_double_R, R_NOMA_EE_double_R, ~, Pi_sys_JT_EE_double_R, Pi_sys_NOMA_EE_double_R, ~, feasible_JT_EE_double_R, feasible_NOMA_EE_double_R, feasible_LOCAL_EE_double_R] = prepare_data_box(data_EE_double_R, aux_kappa, aux_rho, aux_P_fix);
%     [R_JT_EE_double_R_stop_per_BS, R_NOMA_EE_double_R_stop_per_BS, ~, Pi_sys_JT_EE_double_R_stop_per_BS, Pi_sys_NOMA_EE_double_R_stop_per_BS, ~, feasible_JT_EE_double_R_stop_per_BS, feasible_NOMA_EE_double_R_stop_per_BS, feasible_LOCAL_EE_double_R_stop_per_BS] = prepare_data_box(data_EE_double_R_stop_per_BS, aux_kappa, aux_rho, aux_P_fix);
%     [R_JT_EE_half_R_stop_per_BS, R_NOMA_EE_half_R_stop_per_BS, ~, Pi_sys_JT_EE_half_R_stop_per_BS, Pi_sys_NOMA_EE_half_R_stop_per_BS, ~, feasible_JT_EE_half_R_stop_per_BS, feasible_NOMA_EE_half_R_stop_per_BS, feasible_LOCAL_EE_half_R_stop_per_BS] = prepare_data_box(data_EE_half_R_stop_per_BS, aux_kappa, aux_rho, aux_P_fix);
   
%     figure, plot(x_axis,mean(data_EE_half_R.N_iter_local(:,feasible_LOCAL_EE_half_R),2),'-+b','LineWidth',1),
%     hold on 
%     plot(x_axis,mean(data_EE_double_R.N_iter_local(:,feasible_LOCAL_EE_double_R),2),'-*k','LineWidth',1),
%     hold on 
%     plot(x_axis,mean(data_EE_double_R_stop_per_BS.N_iter_local(:,feasible_LOCAL_EE_double_R_stop_per_BS),2),'-*k','LineWidth',1),
%     legend('R/2', '2R','2R stop per BS')
%     xlabel('Minimum data rate requirement (Kbps)');
%     ylabel('Number of iterations');
    non_error_samples_local = ~data_EE_half_R.error_samples.' & ~data_EE_double_R.error_samples.';% & ~data_EE_double_R_stop_per_BS.error_samples.' & ~data_EE_half_R_stop_per_BS.error_samples.';
    feasible_LOCAL_all = feasible_LOCAL_EE_half_R & feasible_LOCAL_EE_double_R & non_error_samples_local;% & feasible_LOCAL_EE_double_R_stop_per_BS & feasible_LOCAL_EE_half_R_stop_per_BS 
    
    feasible_JT = feasible_JT_EE_half_R & feasible_JT_EE_double_R & non_error_samples_local;
    feasible_NOMA = feasible_NOMA_EE_half_R & feasible_NOMA_EE_double_R & non_error_samples_local;
    
    feasible_all = feasible_LOCAL_all & feasible_JT & feasible_NOMA;
    
    N_iter_half = data_EE_half_R.N_iter_local(1,feasible_LOCAL_all);
    N_iter_double = data_EE_double_R.N_iter_local(1,feasible_LOCAL_all);
%     N_iter_double_stop_per_BS = data_EE_double_R_stop_per_BS.N_iter_local(1,feasible_LOCAL_all);
%     N_iter_half_stop_per_BS = data_EE_half_R_stop_per_BS.N_iter_local(1,feasible_LOCAL_all);
    
    figure
    %boxplot([N_iter_half.',N_iter_double.', N_iter_double_stop_per_BS.', N_iter_half_stop_per_BS.'],'Labels',{'Half','Double', 'Double stop per BS', 'Half stop per BS'})
    boxplot([N_iter_half.',N_iter_double.'],'Labels',{'0.5 Mbps','2 Mbps'});
    xlabel('Edge user minimum rate requirement');
    ylabel('Number of iterations');
    saveas(gcf,'fig_review/N_iter.png')
    %save('workspaces/figures_data.mat','Ri_EE','Ri_SE','Ri_min','-append');
    

    ss = 2;
    P_BS1_iter = data_EE_double_R.P_BS1_iter(1,1:50,ss);
    P_BS2_iter = data_EE_double_R.P_BS2_iter(1,1:50,ss);
    EE_total_iter = data_EE_double_R.EE_total_iter(1,1:50,ss);
    
    idxs = find(P_BS1_iter);
    figure, plot(1:length(idxs),P_BS1_iter(idxs),'-+b','LineWidth',1),
    xticks(1:length(idxs))
    xlabel('Iteration index');
    ylabel('Total Power');
    saveas(gcf,'fig_review/powerBS1.png')
    
    figure,
    plot(1:length(idxs),P_BS2_iter(idxs),'-ok','LineWidth',1),
    %legend('BS 1', 'BS 2')
    xticks(1:length(idxs))
    xlabel('Iteration index');
    ylabel('Total Power');
    saveas(gcf,'fig_review/powerBS2.png')
    

    figure, plot(1:length(idxs),EE_total_iter(idxs),'-+b','LineWidth',1),
%     hold on,
%     plot(1:length(idxs),ones(1,length(idxs)).*data_EE_double_R.R_tot_EE_S1_global_1(1,ss)./data_EE_double_R.PC_EE_S1_global(1,ss),'-ok','LineWidth',1),
%      hold on,
%     plot(1:length(idxs),ones(1,length(idxs)).*data_EE_double_R.R_tot_EE_local(1,ss)./data_EE_double_R.PC_EE_local(1,ss),'-*c','LineWidth',1),
    xticks(1:length(idxs))
    xlabel('Iteration index');
    ylabel('EE');
    saveas(gcf,'fig_review/EE_convergence.png')

    % -------------------------------
    % EE - Boxplot
    % -------------------------------
    EE_global_double = data_EE_double_R.R_tot_EE_S1_global_1(1,feasible_LOCAL_all & feasible_JT)./data_EE_double_R.PC_EE_S1_global(1,feasible_LOCAL_all & feasible_JT);
    EE_ILO_double = data_EE_double_R.R_tot_EE_local(1,feasible_LOCAL_all & feasible_JT)./data_EE_double_R.PC_EE_local(1,feasible_LOCAL_all & feasible_JT);
%     figure,
%     boxplot([EE_global.', EE_ILO.'],'Labels',{'Global','ILO'})
    
    EE_global_half = data_EE_half_R.R_tot_EE_S1_global_1(1,feasible_LOCAL_all & feasible_JT)./data_EE_half_R.PC_EE_S1_global(1,feasible_LOCAL_all & feasible_JT);
    EE_ILO_half = data_EE_half_R.R_tot_EE_local(1,feasible_LOCAL_all & feasible_JT)./data_EE_half_R.PC_EE_local(1,feasible_LOCAL_all & feasible_JT);
%     figure,
%     boxplot([EE_global_half.', EE_ILO_half.'],'Labels',{'Global','ILO'})

    EE_global = data_EE.R_tot_EE_S1_global_1(1,feasible_LOCAL_all & feasible_JT)./data_EE.PC_EE_S1_global(1,feasible_LOCAL_all & feasible_JT);
    EE_ILO = data_EE.R_tot_EE_local(1,feasible_LOCAL_all & feasible_JT)./data_EE.PC_EE_local(1,feasible_LOCAL_all & feasible_JT);

    box_data = {[EE_global_half; EE_global_double].', [EE_ILO_half; EE_ILO_double].'};
    figure,boxplotGroup(box_data, 'PrimaryLabels', {'Global JTCN' 'ILO'}, ...
      'SecondaryLabels',cellstr(["Edge user R^{min}=0.5 Mbps","Edge user R^{min}=2 Mbps"]), 'InterGroupSpace', 1, ...
       'GroupType','withinGroups')
    ylabel('Network EE')
    %xlabel('Edge user minimum rate requirement');
    saveas(gcf,'fig_review/EE_global_vs_ILO.png')
    mean(EE_global_half)./mean(EE_ILO_half)
    mean(EE_global_double)./mean(EE_ILO_double)
    mean(EE_global)./mean(EE_ILO)
    
%     EE_global_JT = data_EE_double_R.R_tot_EE_S1_global_1(1,feasible_all)./data_EE_double_R.PC_EE_S1_global(1,feasible_all);
%     EE_global_NOMA = data_EE_double_R.R_tot_EE_S3_global_1(1,feasible_all)./data_EE_double_R.PC_EE_S3_global(1,feasible_all);
%     EE_ILO = data_EE_double_R.R_tot_EE_local(1,feasible_all)./data_EE_double_R.PC_EE_local(1,feasible_all);
%     figure,
%     boxplot([EE_global_JT.', EE_ILO.', EE_global_NOMA.'],'Labels',{'Global JT','ILO', 'Global NOMA'})

    EE_JTCN_half = data_EE_half_R.R_tot_EE_S1_global_1(1,feasible_all)./data_EE_half_R.PC_EE_S1_global(1,feasible_all);
    EE_JTCN_double = data_EE_double_R.R_tot_EE_S1_global_1(1,feasible_all)./data_EE_double_R.PC_EE_S1_global(1,feasible_all);
    EE_NOMA_half = data_EE_half_R.R_tot_EE_S3_global_1(1,feasible_all)./data_EE_half_R.PC_EE_S3_global(1,feasible_all);
    EE_NOMA_double = data_EE_double_R.R_tot_EE_S3_global_1(1,feasible_all)./data_EE_double_R.PC_EE_S3_global(1,feasible_all);
    EE_ILO_half = data_EE_half_R.R_tot_EE_local(1,feasible_all)./data_EE_half_R.PC_EE_local(1,feasible_all);
    EE_ILO_double = data_EE_double_R.R_tot_EE_local(1,feasible_all)./data_EE_double_R.PC_EE_local(1,feasible_all);
    
box_data = {[EE_JTCN_half; EE_JTCN_double].', [EE_NOMA_half; EE_NOMA_double].', [EE_ILO_half; EE_ILO_double].'};
    figure,boxplotGroup(box_data, 'PrimaryLabels', {'JTCN' 'NOMA','ILO'}, ...
      'SecondaryLabels',cellstr(["Edge user R^{min}=0.5 Mbps","Edge user R^{min}=2 Mbps"]), 'InterGroupSpace', 1, ...
       'GroupType','withinGroups')
    ylabel('Network EE')
    saveas(gcf,'fig_review/EE_JTCN_vs_NOMA.png')
    mean(EE_JTCN_half)./mean(EE_NOMA_half)
    mean(EE_JTCN_double)./mean(EE_NOMA_double)
    
    mean(EE_JTCN_half)./mean(EE_ILO_half)
    mean(EE_JTCN_double)./mean(EE_ILO_double)
    
    save('workspaces/new_figs_data.mat','EE_JTCN_half','EE_JTCN_double','EE_NOMA_half','EE_NOMA_double','EE_ILO_half','EE_ILO_double', 'N_iter_half', 'N_iter_double');
   
    
    
    
 