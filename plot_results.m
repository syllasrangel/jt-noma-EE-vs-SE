function [] = plot_results(filepath, plot_outage, plot_iter, plot_throughput, plot_EE, plot_EC, plot_individual_EE)
% plot_outage = false;
% plot_iter = false;
% plot_throughput = false;
% plot_EE = false;
% plot_EC = false;
% plot_individual_EE = true;

data = load(filepath);

data.kappa_aux = 0.5; % remove this.

%data_kappa_250 = load("workspaces\2022_11_19_14_04_workspace_PCM_proposed_kappa_2_50_rho_0_10_100samp.mat");
%data_kappa_050 =
%load("workspaces\2022_11_19_14_04_workspace_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat"); %With local
%data_kappa_050 = load("workspaces\2022_11_27_23_00_workspace_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat"); %without local
%data_kappa_00 = load("workspaces\2022_11_19_14_03_workspace_PCM_M2_kappa_0_00_rho_0_10_100samp.mat");



%data_PCM1 = load("workspaces\2022_11_23_10_11_workspace_PCM_M1_kappa_0_00_rho_0_00_100samp.mat");
%data_PCM2 = load("workspaces\2022_11_26_00_30_workspace_PCM_M2_kappa_0_00_rho_0_00_100samp.mat");
data_PCM1 = load("workspaces\2022_12_08_18_31_workspace_PCM_M1_kappa_0_00_rho_0_00_100samp.mat");
data_PCM2 = load("workspaces\2022_12_09_00_25_workspace_PCM_M2_kappa_0_00_rho_0_00_100samp.mat");
data_kappa_00 = load("workspaces\2022_12_12_15_50_workspace_PCM_proposed_kappa_0_00_rho_0_10_100samp.mat"); % Without local
data_kappa_050 = load('workspaces\2022_12_08_18_32_workspace_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat');
%data_kappa_250 = load("workspaces\2022_12_11_23_16_workspace_PCM_proposed_kappa_2_50_rho_0_10_100samp.mat"); % Without local
data_kappa_250 = load("workspaces\2022_12_11_23_18_workspace_PCM_proposed_kappa_2_50_rho_0_10_100samp.mat");




% w=5
% data_PCM1 = load("workspaces\2022_12_01_00_00_workspace_PCM_M1_kappa_0_00_rho_0_00_100samp.mat");
% data_PCM2 = load("workspaces\2022_12_01_00_42_workspace_PCM_M2_kappa_0_00_rho_0_00_100samp.mat");
% data_kappa_050 = load("workspaces\2022_11_30_23_49_workspace_PCM_proposed_kappa_0_50_rho_0_10_100samp.mat"); 

x_axis = data.x_axis;
N_samples = data.N_samples;

if (data.PCM == "M1")
    PCM_text = "PCM 1";
elseif (data.PCM == "M2")
    PCM_text = "PCM 2";
elseif  (data.PCM == "proposed")
    PCM_text = "PCM proposed";
end


if(data.s~=N_samples)
    interval = 1:data.s-1;
else
    interval = 1:data.s;
end
simu_samps = zeros(1,N_samples);
simu_samps(interval) = 1;
% & sum(isnan(data.Exit_EE_S3_global),1)==0 excludes the samples with
% NaN in the exit flag (only happens in the skiped samples)
feasible_JT = (~sum((data.Exit_EE_S1_global < 0),1)>0) & simu_samps & sum(isnan(data.Exit_EE_S1_global),1)==0;
feasible_NOMA = (~sum((data.Exit_EE_S3_global < 0),1)>0) & simu_samps & sum(isnan(data.Exit_EE_S3_global),1)==0;
feasible_LOCAL = (~sum((data.Exit_EE_local < 0),1)>0) & simu_samps & sum(isnan(data.Exit_EE_local),1)==0;

% Remove the outage scenarios from all simulations
feasible_all = feasible_JT & feasible_NOMA & feasible_LOCAL;


if(plot_outage)
    figure, 
    plot(x_axis,mean(data.Exit_EE_S3_global(:,simu_samps & sum(isnan(data.Exit_EE_S3_global),1)==0) < 0,2),'-.ok','LineWidth',1),
    hold on,
    plot(x_axis,mean(data.Exit_EE_S1_global(:,simu_samps & sum(isnan(data.Exit_EE_S1_global),1)==0) < 0,2),'-+b','LineWidth',1),
    hold on,
    plot(x_axis,mean(data.Exit_EE_local(:,simu_samps & sum(isnan(data.Exit_EE_local),1)==0) < 0,2),'--og','LineWidth',1)    
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Outage rate');
    title('Outage rate');
    legend('Conventional NOMA', 'JT-CoMP NOMA', 'Local DPS-CoMP NOMA','Location', 'southeast')
end

% overall_outage = (Exit_EE_S1 <= 0) + (Exit_EE_S3 <= 0) + (Exit_S1 <= 0) + (Exit_S3 <= 0) > 0;
% figure, plot(x_axis,mean(overall_outage,2),'-+b','LineWidth',1)
% xlabel('Minimum data rate requirement (Kbps)');
% ylabel('Outage probability');
% title('Overall Outage');


if (plot_iter)
    figure, plot(x_axis,mean(data.N_iter_S3_global(:,feasible_all),2),'-.ok','LineWidth',1),
    hold on,
    plot(x_axis,mean(data.N_iter_S1_global(:,feasible_all),2),'-+b','LineWidth',1),
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Number of iterations');
    legend('Conventional NOMA', 'JT-CoMP NOMA','Location', 'northwest')
    title(sprintf('Dinkelbach - %s', PCM_text));
    
    figure, plot(x_axis,mean(data.N_iter_local(:,feasible_all),2),'-.ok','LineWidth',1),
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Number of iterations');
end

R_EE_S1_global = data.R_EE_S1_global;
R_EE_S3_global = data.R_EE_S3_global;
R_EE_local = data.R_EE_local;
R_EE_S1_global(isnan(R_EE_S1_global)) = 0;
R_EE_S3_global(isnan(R_EE_S3_global)) = 0;
R_EE_local(isnan(R_EE_local)) = 0;
R_tot_EE_S1_global = squeeze(sum(sum(R_EE_S1_global,2),3));
R_tot_EE_S3_global = squeeze(sum(sum(R_EE_S3_global,2),3));
R_tot_EE_local = squeeze(sum(sum(R_EE_local,2),3));



if(plot_throughput)
    [R_EE_t_S1_global, R_EE_t_CI_S1_global] = mean_confidence_interval(R_tot_EE_S1_global(:,feasible_all));
    [R_EE_t_S3_global, R_EE_t_CI_S3_global] = mean_confidence_interval(R_tot_EE_S3_global(:,feasible_all));
    [R_EE_t_local, R_EE_t_CI_local] = mean_confidence_interval(R_tot_EE_local(:,feasible_all));
    
    % ALL
    figure;
    hAx=axes;
    errorbar(x_axis,R_EE_t_S1_global, R_EE_t_CI_S1_global(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,R_EE_t_S3_global,R_EE_t_CI_S3_global(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,R_EE_t_local,R_EE_t_CI_local(:,2),'-.og','LineWidth',1),
    hold on,
    legend('EE JT-CoMP-NOMA','EE Conventional NOMA','ILO');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Average throughput (b/s)');
    title(sprintf('Throughput - %s',PCM_text));
    hAx.YScale='log';
    
    %JT only
    [R_EE_t_S1_global, R_EE_t_CI_S1_global] = mean_confidence_interval(R_tot_EE_S1_global(:,feasible_JT));
    figure;
    hAx=axes;
    errorbar(x_axis,R_EE_t_S1_global, R_EE_t_CI_S1_global(:,2),'-+b','LineWidth',1)
    legend('EE JT-CoMP-NOMA');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Average throughput (b/s)');
    title(sprintf('Throughput - %s',PCM_text));
    hAx.YScale='log';
end

% %figure, boxplot(R_tot_EE_S3_global(:,samples).', x_axis)
% 
% bpColors = cool(length(x_axis));
% %R_tot_EE_S1_global(R_tot_EE_S1_global == 0) = NaN;
% %R_tot_EE_S3_global(R_tot_EE_S3_global == 0) = NaN;
% data = {R_tot_EE_S1_global(:,samples).', R_tot_EE_S3_global(:,samples).'};
% figure,boxplotGroup(data, 'PrimaryLabels', {'JT' 'Conv.'}, ...
%   'SecondaryLabels',cellstr(string(x_axis)), 'InterGroupSpace', 1, ...
%   'BoxStyle','filled',...
%     'Colors',bpColors,'GroupType','withinGroups')
% title('Global solution - Average throughput');
% 
% 
% 
% % OBS.: It seems that forcing the edge user to be in the first BS cluster 
% % makes the samples with both JT and Conv feasible be only the ones where
% % the channel BS1 -> edge user is good.
% figure,
% hold on,
% xlabel('gamma BS 1 to edge user');
% ylabel('gamma BS 2 to edge user');
% title('Feasibility for all R_{min} for both JT and conv.')
% legended1 = false;
% legended2 = false;
% for samp = 1:N_samples
%     if(samples(samp))
%         if(legended1 == true)
%             h = plot(gamma_values(N_inner_users+N_JT_users,1,1,samp),gamma_values(N_inner_users+N_JT_users,2,2,samp),'g.');
%             h.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         else
%             plot(gamma_values(N_inner_users+N_JT_users,1,1,samp),gamma_values(N_inner_users+N_JT_users,2,2,samp),'g.', 'DisplayName','Feasible for all R_{min}')
%         end
%         legended1 = true;
%     else
%         if(legended2 == true)
%             h = plot(gamma_values(N_inner_users+N_JT_users,1,1,samp),gamma_values(N_inner_users+N_JT_users,2,2,samp),'r.');
%             h.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         else
%             plot(gamma_values(N_inner_users+N_JT_users,1,1,samp),gamma_values(N_inner_users+N_JT_users,2,2,samp),'r.', 'DisplayName','Infeasible for some R_{min}')
%         end
%         legended2 = true;
%     end
% end
% hold off
% lgd = legend;
% 
% 
% % Feasibility for JT
% figure,
% hold on,
% xlabel('gamma BS 1 to edge user');
% ylabel('gamma BS 2 to edge user');
% title('Feasibility for all R_{min} for JT')
% legended1 = false;
% legended2 = false;
% for samp = 1:N_samples
%     if(sum(Exit_EE_S1_global(:,samp)~=1)==0)
%         if(legended1 == true)
%             h = plot(gamma_values(N_inner_users+N_JT_users,1,1,samp),gamma_values(N_inner_users+N_JT_users,2,2,samp),'g.');
%             h.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         else
%             plot(gamma_values(N_inner_users+N_JT_users,1,1,samp),gamma_values(N_inner_users+N_JT_users,2,2,samp),'g.', 'DisplayName','Feasible for all R_{min}')
%         end
%         legended1 = true;
%     else
%         if(legended2 == true)
%             h = plot(gamma_values(N_inner_users+N_JT_users,1,1,samp),gamma_values(N_inner_users+N_JT_users,2,2,samp),'r.');
%             h.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         else
%             plot(gamma_values(N_inner_users+N_JT_users,1,1,samp),gamma_values(N_inner_users+N_JT_users,2,2,samp),'r.', 'DisplayName','Infeasible for some R_{min}')
%         end
%         legended2 = true;
%     end
% end
% hold off
% lgd = legend;
% 
% % Feasibility for Conventional
% figure,
% hold on,
% xlabel('gamma BS 1 to edge user');
% ylabel('gamma BS 2 to edge user');
% title('Feasibility for all R_{min} for conv.')
% legended1 = false;
% legended2 = false;
% legended3 = false;
% for samp = 1:N_samples
%     if(sum(Exit_EE_S3_global(:,samp)~=1)==0)
%         if(Pi_EE_S1_global(1,3,samp)==0)
%             if(sum(Pi_EE_S1_global(:,3,samp))~=0) % Both BSs transmit with power ~= 0
%                 color = 'y.';
%             else
%                 color = 'b.';
%             end
%         else
%             color = 'g.';
%         end
%         if(legended1 == true)
%             if(Pi_EE_S1_global(1,3,samp)~=0 || legended3)
%                 h = plot(gamma_values(N_inner_users+N_JT_users,1,1,samp),gamma_values(N_inner_users+N_JT_users,2,2,samp),color);
%                 h.Annotation.LegendInformation.IconDisplayStyle = 'off';
%             else
%                 plot(gamma_values(N_inner_users+N_JT_users,1,1,samp),gamma_values(N_inner_users+N_JT_users,2,2,samp),color, 'DisplayName','Feasible for all R_{min} | P = 0 at BS1 for CoMP user')
%                 legended3 = true;
%             end
%         else
%                 plot(gamma_values(N_inner_users+N_JT_users,1,1,samp),gamma_values(N_inner_users+N_JT_users,2,2,samp),color, 'DisplayName','Feasible for all R_{min}')
%                 legended1 = true;
%         end
%         
%     else
%         if(legended2 == true)
%             h = plot(gamma_values(N_inner_users+N_JT_users,1,1,samp),gamma_values(N_inner_users+N_JT_users,2,2,samp),'r.');
%             h.Annotation.LegendInformation.IconDisplayStyle = 'off';
%         else
%             plot(gamma_values(N_inner_users+N_JT_users,1,1,samp),gamma_values(N_inner_users+N_JT_users,2,2,samp),'r.', 'DisplayName','Infeasible for some R_{min}')
%         end
%         legended2 = true;
%     end
% end
% hold off
% lgd = legend;

% % power difference S1 vs S3
% samples2 = (samples - (squeeze(Pi_EE_S1_global(1,3,:)==0).' & samples))==1;
% P_diff = (abs(Pi_EE_S3_global(:,:,samples2) - Pi_EE_S1_global(:,1:5,samples2))./Pi_EE_S1_global(:,1:5,samples2))*100;
% 
% figure,
% plot(x_axis, mean(mean(P_diff,3),2))
% xlabel('Minimum data rate requirement (Kbps)');
% ylabel('Average power difference (%)');
% title('Difference in the allocated power for JT and conv.')


% Computes the system power expenditure
Pi_sys_EE_S1_global_PCM1 = zeros(length(x_axis),N_samples);
Pi_sys_EE_S3_global_PCM1 = zeros(length(x_axis),N_samples);
Pi_sys_EE_local_PCM1 = zeros(length(x_axis),N_samples);
Pi_sys_EE_S1_global_prop = zeros(length(x_axis),N_samples);
Pi_sys_EE_S3_global_prop = zeros(length(x_axis),N_samples);
Pi_sys_EE_local_prop = zeros(length(x_axis),N_samples);
Pi_sys_EE_S1_global_PCM2 = zeros(length(x_axis),N_samples);
Pi_sys_EE_S3_global_PCM2 = zeros(length(x_axis),N_samples);
Pi_sys_EE_local_PCM2 = zeros(length(x_axis),N_samples);
for r=1:length(data.R_Kbps)
    for ss=1:N_samples
        Pi_sys_EE_S1_global_PCM1(r,ss) = system_power_consumption(data.Pi_EE_S1_global(r,:,ss), data.gamma, 0, 0, 0, true, false);
        Pi_sys_EE_S3_global_PCM1(r,ss) = system_power_consumption(data.Pi_EE_S3_global(r,:,ss), data.gamma, 0, 0, 0, false, false);
        Pi_sys_EE_local_PCM1(r,ss) = system_power_consumption(data.Pi_EE_local(r,:,ss), data.gamma, 0, 0, 0, false, false);
        
        % TODO: Do I consider rho here? Check in the paper
%         Pi_sys_EE_S1_global_PCM2(r,ss) = system_power_consumption(data.Pi_EE_S1_global(r,:,ss), data.gamma, data.rho_aux, data.P_fix_aux, 0, true, false);
%         Pi_sys_EE_S3_global_PCM2(r,ss) = system_power_consumption(data.Pi_EE_S3_global(r,:,ss), data.gamma, data.rho_aux, data.P_fix_aux, 0, false, false);
%         Pi_sys_EE_local_PCM2(r,ss) = system_power_consumption(data.Pi_EE_local(r,:,ss), data.gamma, data.rho_aux, data.P_fix_aux, 0, false, false);
        Pi_sys_EE_S1_global_PCM2(r,ss) = system_power_consumption(data.Pi_EE_S1_global(r,:,ss), data.gamma, 0, data.P_fix_aux, 0, true, false);
        Pi_sys_EE_S3_global_PCM2(r,ss) = system_power_consumption(data.Pi_EE_S3_global(r,:,ss), data.gamma, 0, data.P_fix_aux, 0, false, false);
        Pi_sys_EE_local_PCM2(r,ss) = system_power_consumption(data.Pi_EE_local(r,:,ss), data.gamma, 0, data.P_fix_aux, 0, false, false);
        
        Pi_sys_EE_S1_global_prop(r,ss) = system_power_consumption(data.Pi_EE_S1_global(r,:,ss), data.gamma, data.rho_aux, data.P_fix_aux, data.kappa_aux, true, false);
        Pi_sys_EE_S3_global_prop(r,ss) = system_power_consumption(data.Pi_EE_S3_global(r,:,ss), data.gamma, data.rho_aux, data.P_fix_aux, data.kappa_aux, false, false);
        Pi_sys_EE_local_prop(r,ss) = system_power_consumption(data.Pi_EE_local(r,:,ss), data.gamma, data.rho_aux, data.P_fix_aux, data.kappa_aux, false, false);
    end
end 



EE_var_kappa = true;
if(EE_var_kappa)    
    %figure, boxplot(R_tot_EE_S3_global(:,feasible_all).', x_axis)
    
    [R_JT_00, R_NOMA_00, R_LOCAL_00, Pi_sys_JT_00, Pi_sys_NOMA_00, Pi_sys_local_00, feasible_JT_00, feasible_NOMA_00, feasible_LOCAL_00] = prepare_data_box(data_kappa_00);
    [R_JT_050, R_NOMA_050, R_LOCAL_050, Pi_sys_JT_050, Pi_sys_NOMA_050, Pi_sys_local_050, feasible_JT_050, feasible_NOMA_050, feasible_LOCAL_050] = prepare_data_box(data_kappa_050);
    [R_JT_250, R_NOMA_250, R_LOCAL_250, Pi_sys_JT_250, Pi_sys_NOMA_250, Pi_sys_local_250, feasible_JT_250, feasible_NOMA_250, feasible_LOCAL_250] = prepare_data_box(data_kappa_250);
    
    R_idx = 4;
    feasible = feasible_JT_00 & feasible_JT_050 & feasible_JT_250 & feasible_NOMA_00 & feasible_NOMA_050 & feasible_NOMA_250;
    EE_00 = (R_JT_00(R_idx,feasible)./Pi_sys_JT_00(R_idx,feasible)).';
    EE_050 = (R_JT_050(R_idx,feasible)./Pi_sys_JT_050(R_idx,feasible)).';
    EE_250 = (R_JT_250(R_idx,feasible)./Pi_sys_JT_250(R_idx,feasible)).';
    NOMA_00 = (R_NOMA_00(R_idx,feasible)./Pi_sys_NOMA_00(R_idx,feasible)).';
    NOMA_050 = (R_NOMA_050(R_idx,feasible)./Pi_sys_NOMA_050(R_idx,feasible)).';
    NOMA_250 = (R_NOMA_250(R_idx,feasible)./Pi_sys_NOMA_250(R_idx,feasible)).';
    %min_size = min([length(EE_00), length(EE_050), length(EE_250)]);
    
    %figure, boxplot([EE_00(1:min_size), EE_050(1:min_size), EE_250(1:min_size)], ["kappa = 0","kappa = 0.5","kappa = 2.5"])
    
    box_data = {[EE_00, EE_050, EE_250], [NOMA_00, NOMA_050, NOMA_250]};
    figure,boxplotGroup(box_data, 'PrimaryLabels', {'JT' 'Conv.'}, ...
      'SecondaryLabels',cellstr(["\kappa = 0","\kappa = 0.5","\kappa = 2.5"]), 'InterGroupSpace', 1, ...
       'GroupType','withinGroups')
    title('Energy Efficiency');
    ylabel('System energy efficiency (b/s/Joule)')
    
    
    
%     R_idx = 4;
%     EE_00 = (R_JT_00(R_idx,feasible_JT_00)./Pi_sys_JT_00(R_idx,feasible_JT_00)).';
%     EE_050 = (R_JT_050(R_idx,feasible_JT_050)./Pi_sys_JT_050(R_idx,feasible_JT_050)).';
%     EE_250 = (R_JT_250(R_idx,feasible_JT_250)./Pi_sys_JT_250(R_idx,feasible_JT_250)).';
%     NOMA_00 = (R_NOMA_00(R_idx,feasible_NOMA_00)./Pi_sys_NOMA_00(R_idx,feasible_NOMA_00)).';
%     NOMA_050 = (R_NOMA_050(R_idx,feasible_NOMA_050)./Pi_sys_NOMA_050(R_idx,feasible_NOMA_050)).';
%     NOMA_250 = (R_NOMA_250(R_idx,feasible_NOMA_250)./Pi_sys_NOMA_250(R_idx,feasible_NOMA_250)).';
%     min_size = min([length(EE_00), length(EE_050), length(EE_250), length(NOMA_00), length(NOMA_050), length(NOMA_250)]);
%     
%     box_data = {[EE_00(1:min_size), EE_050(1:min_size), EE_250(1:min_size)], [NOMA_00(1:min_size), NOMA_050(1:min_size), NOMA_250(1:min_size)]};
%     figure,boxplotGroup(box_data, 'PrimaryLabels', {'JT' 'Conv.'}, ...
%       'SecondaryLabels',cellstr(["\kappa = 0","\kappa = 0.5","\kappa = 2.5"]), 'InterGroupSpace', 1, ...
%        'GroupType','withinGroups')
%     title('Energy Efficiency');
%     ylabel('System energy efficiency (b/s/Joule)')
end

EE_PCMs = true;
if(EE_PCMs)    
    aux_kappa = 0.5;
    aux_rho = 0.1;
    aux_P_fix = 1;
    [R_JT_PCM1, R_NOMA_PCM1, R_LOCAL_PCM1, Pi_sys_JT_PCM1, Pi_sys_NOMA_PCM1, Pi_sys_local_PCM1, feasible_JT_PCM1, feasible_NOMA_PCM1, feasible_LOCAL_PCM1] = prepare_data_box(data_PCM1, aux_kappa, aux_rho, aux_P_fix);
    [R_JT_PCM2, R_NOMA_PCM2, R_LOCAL_PCM2, Pi_sys_JT_PCM2, Pi_sys_NOMA_PCM2, Pi_sys_local_PCM2, feasible_JT_PCM2, feasible_NOMA_PCM2, feasible_LOCAL_PCM2] = prepare_data_box(data_PCM2, aux_kappa, aux_rho, aux_P_fix);
    [R_JT_050, R_NOMA_050, R_LOCAL_050, Pi_sys_JT_050, Pi_sys_NOMA_050, Pi_sys_local_050, feasible_JT_050, feasible_NOMA_050, feasible_LOCAL_050] = prepare_data_box(data_kappa_050, aux_kappa, aux_rho, aux_P_fix);
    %[R_JT_250, R_NOMA_250, R_LOCAL_250, Pi_sys_JT_250, Pi_sys_NOMA_250, Pi_sys_local_250, feasible_JT_250, feasible_NOMA_250, feasible_LOCAL_250] = prepare_data_box(data_kappa_250, aux_kappa, aux_rho, aux_P_fix);
    
    feasible = feasible_JT_PCM1 & feasible_JT_PCM2 & feasible_JT_050; %& feasible_JT_250;
    [EE_EE_t_JT_PCM1, EE_EE_t_CI_JT_PCM1] = mean_confidence_interval(R_JT_PCM1(:,feasible)./Pi_sys_JT_PCM1(:,feasible));
    [EE_EE_t_JT_PCM2, EE_EE_t_CI_JT_PCM2] = mean_confidence_interval(R_JT_PCM2(:,feasible)./Pi_sys_JT_PCM2(:,feasible));
    [EE_EE_t_JT_prop, EE_EE_t_CI_JT_prop] = mean_confidence_interval(R_JT_050(:,feasible)./Pi_sys_JT_050(:,feasible));
    %[EE_EE_t_JT_prop_250, EE_EE_t_CI_JT_prop_250] = mean_confidence_interval(R_JT_250(:,feasible)./Pi_sys_JT_250(:,feasible));
    
%     feasible = feasible_NOMA_PCM1 & feasible_NOMA_PCM2 & feasible_NOMA_050;
%     [EE_EE_t_JT_PCM1, EE_EE_t_CI_JT_PCM1] = mean_confidence_interval(R_NOMA_PCM1(:,feasible)./Pi_sys_NOMA_PCM1(:,feasible));
%     [EE_EE_t_JT_PCM2, EE_EE_t_CI_JT_PCM2] = mean_confidence_interval(R_NOMA_PCM2(:,feasible)./Pi_sys_NOMA_PCM2(:,feasible));
%     [EE_EE_t_JT_prop, EE_EE_t_CI_JT_prop] = mean_confidence_interval(R_NOMA_050(:,feasible)./Pi_sys_NOMA_050(:,feasible));
    
    figure,
    hAx=axes;
    errorbar(x_axis,EE_EE_t_JT_PCM1, EE_EE_t_CI_JT_PCM1(:,2),'-+b','LineWidth',1),
    hold on,
    errorbar(x_axis,EE_EE_t_JT_PCM2, EE_EE_t_CI_JT_PCM2(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,EE_EE_t_JT_prop, EE_EE_t_CI_JT_prop(:,2),'-.og','LineWidth',1),
    %hold on,
    %errorbar(x_axis,EE_EE_t_JT_prop_250, EE_EE_t_CI_JT_prop_250(:,2),'-.og','LineWidth',1),
    hold on,
    %legend('Opt. PCM 1','Opt. PCM 2','Opt. PCM prop. \kappa = 0.50','Opt. PCM prop. \kappa = 2.50');
    legend('Opt. PCM 1','Opt. PCM 2','Opt. PCM prop. \kappa = 0.50');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Average energy efficiency (b/s/Joule)');
    title(sprintf('JT - Evaluated with PCM prop (rho = %0.2f, kappa = %0.2f, P_{fix} = %0.2f)', aux_rho, aux_kappa, aux_P_fix));
    %hAx.YScale='log';
    
    [TP_EE_t_JT_PCM1, TP_EE_t_CI_JT_PCM1] = mean_confidence_interval(R_JT_PCM1(:,feasible));
    [TP_EE_t_JT_PCM2, TP_EE_t_CI_JT_PCM2] = mean_confidence_interval(R_JT_PCM2(:,feasible));
    [TP_EE_t_JT_prop, TP_EE_t_CI_JT_prop] = mean_confidence_interval(R_JT_050(:,feasible));
    %[EE_EE_t_JT_prop_250, EE_EE_t_CI_JT_prop_250] = mean_confidence_interval(R_JT_250(:,feasible));
    
    figure,
    hAx=axes;
    errorbar(x_axis,TP_EE_t_JT_PCM1, TP_EE_t_CI_JT_PCM1(:,2),'-+b','LineWidth',1),
    hold on,
    errorbar(x_axis,TP_EE_t_JT_PCM2, TP_EE_t_CI_JT_PCM2(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,TP_EE_t_JT_prop, TP_EE_t_CI_JT_prop(:,2),'-.og','LineWidth',1),
    hold on,
    plot(x_axis, (data_PCM1.N_BSs*data_PCM1.N_inner_users + data_PCM1.N_JT_users)*data_PCM1.R,'-r','LineWidth',1)
    legend('Opt. PCM 1','Opt. PCM 2','Opt. PCM prop. \kappa = 0.50');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Average throughput (b/s)');
    title(sprintf('JT - Evaluated with PCM prop (rho = %0.2f, kappa = %0.2f, P_{fix} = %0.2f)', aux_rho, aux_kappa, aux_P_fix));
    %hAx.YScale='log';
    
    
    aloc_power_JT_PCM1 = squeeze(sum(data_PCM1.Pi_EE_S1_global,2));
    aloc_power_NOMA_PCM1 = squeeze(sum(data_PCM1.Pi_EE_S3_global,2));
    aloc_power_LOCAL_PCM1 = squeeze(sum(data_PCM1.Pi_EE_local,2));

    aloc_power_JT_PCM2 = squeeze(sum(data_PCM2.Pi_EE_S1_global,2));
    aloc_power_NOMA_PCM2 = squeeze(sum(data_PCM2.Pi_EE_S3_global,2));
    aloc_power_LOCAL_PCM2 = squeeze(sum(data_PCM2.Pi_EE_local,2));
    
    aloc_power_JT_PCMprop = squeeze(sum(data_kappa_050.Pi_EE_S1_global,2));
    aloc_power_NOMA_PCMprop = squeeze(sum(data_kappa_050.Pi_EE_S3_global,2));
    aloc_power_LOCAL_PCMprop = squeeze(sum(data_kappa_050.Pi_EE_local,2));
            
              
    [avg_aloc_power_JT_PCM1, ~] = mean_confidence_interval(aloc_power_JT_PCM1(:, feasible_JT_PCM1));
    [avg_aloc_power_JT_PCM2, ~] = mean_confidence_interval(aloc_power_JT_PCM2(:, feasible_JT_PCM2));
    [avg_aloc_power_JT_PCMprop, ~] = mean_confidence_interval(aloc_power_JT_PCMprop(:, feasible_JT_050));
    
    
    rates = [1,4,7];
    box_data = {(R_tot_EE_S1_global(rates,feasible_all)./Pi_sys_EE_S1_global_prop(rates,feasible_all)).', (R_tot_EE_S3_global(rates,feasible_all)./Pi_sys_EE_S3_global_prop(rates,feasible_all)).', (R_tot_EE_local(rates,feasible_all)./Pi_sys_EE_local_prop(rates,feasible_all)).'};
    figure,boxplotGroup(box_data, 'PrimaryLabels', {'JT' 'Conv.', 'ILO'}, ...
      'SecondaryLabels',cellstr(string(x_axis(rates))), 'InterGroupSpace', 1, ...
       'GroupType','withinGroups')
    title('Energy Efficiency');
    ylabel('System energy efficiency (b/s/Joule)')
    
    
    rates = [1,4,7];
    box_data = {(R_tot_EE_S1_global(rates,feasible_all)).', (R_tot_EE_S3_global(rates,feasible_all)).', (R_tot_EE_local(rates,feasible_all)).'};
    figure,boxplotGroup(box_data, 'PrimaryLabels', {'JT' 'Conv.', 'ILO'}, ...
      'SecondaryLabels',cellstr(string(x_axis(rates))), 'InterGroupSpace', 1, ...
       'GroupType','withinGroups')
    %title('Throughput');
    ylabel('Average throughput (b/s)')
    
end




if(plot_EE)
    
    [EE_EE_t_S1_global_PCM_prop, EE_EE_t_CI_S1_global_PCM_prop] = mean_confidence_interval(R_tot_EE_S1_global(:,feasible_all)./Pi_sys_EE_S1_global_prop(:,feasible_all));
    [EE_EE_t_S3_global_PCM_prop, EE_EE_t_CI_S3_global_PCM_prop] = mean_confidence_interval(R_tot_EE_S3_global(:,feasible_all)./Pi_sys_EE_S3_global_prop(:,feasible_all));
    [EE_EE_t_local_PCM_prop, EE_EE_t_CI_local_PCM_prop] = mean_confidence_interval(data.R_tot_EE_local(:,feasible_all)./Pi_sys_EE_local_prop(:,feasible_all));

    figure;
    hAx=axes;
    errorbar(x_axis,EE_EE_t_S1_global_PCM_prop, EE_EE_t_CI_S1_global_PCM_prop(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,EE_EE_t_S3_global_PCM_prop, EE_EE_t_CI_S3_global_PCM_prop(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,EE_EE_t_local_PCM_prop, EE_EE_t_CI_local_PCM_prop(:,2),'-.og','LineWidth',1),
    hold on,
    legend('EE JT-CoMP-NOMA','EE Conventional NOMA','ILO');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Average energy efficiency (b/s/Joule)');
    title(sprintf('EE - Optimized with %s (rho = %0.2f, kappa = %0.2f, P_{fix} = %0.2f) evaluated with PCM prop (rho = %0.2f, kappa = %0.2f, P_{fix} = %0.2f)',PCM_text, data.rho, data.kappa, data.P_fix,data.rho_aux, data.kappa_aux, data.P_fix_aux));
    hAx.YScale='log';
    
    [EE_EE_t_S1_global, EE_EE_t_CI_S1_global] = mean_confidence_interval(R_tot_EE_S1_global(:,feasible_all)./Pi_sys_EE_S1_global_PCM2(:,feasible_all));
    [EE_EE_t_S3_global, EE_EE_t_CI_S3_global] = mean_confidence_interval(R_tot_EE_S3_global(:,feasible_all)./Pi_sys_EE_S3_global_PCM2(:,feasible_all));
    [EE_EE_t_local, EE_EE_t_CI_local] = mean_confidence_interval(R_tot_EE_local(:,feasible_all)./Pi_sys_EE_local_PCM2(:,feasible_all));

    figure;
    hAx=axes;
    errorbar(x_axis,EE_EE_t_S1_global, EE_EE_t_CI_S1_global(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,EE_EE_t_S3_global, EE_EE_t_CI_S3_global(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,EE_EE_t_local, EE_EE_t_CI_local(:,2),'-.og','LineWidth',1),
    hold on,
    legend('EE JT-CoMP-NOMA','EE Conventional NOMA','ILO');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Average energy efficiency (b/s/Joule)');
    title(sprintf('EE - Optimized with %s (rho = %0.2f, kappa = %0.2f, P_{fix} = %0.2f) evaluated with PCM 2 (rho = %0.2f, kappa = %0.2f, P_{fix} = %0.2f)',PCM_text, data.rho, data.kappa, data.P_fix, 0, 0, data.P_fix_aux));
    hAx.YScale='log';
    
%     if(data.PCM ~= "proposed")
    [EE_EE_t_S1_global, EE_EE_t_CI_S1_global] = mean_confidence_interval(R_tot_EE_S1_global(:,feasible_all)./Pi_sys_EE_S1_global_PCM1(:,feasible_all));
    [EE_EE_t_S3_global, EE_EE_t_CI_S3_global] = mean_confidence_interval(R_tot_EE_S3_global(:,feasible_all)./Pi_sys_EE_S3_global_PCM1(:,feasible_all));
    [EE_EE_t_local, EE_EE_t_CI_local] = mean_confidence_interval(R_tot_EE_local(:,feasible_all)./Pi_sys_EE_local_PCM1(:,feasible_all));

    figure;
    hAx=axes;
    errorbar(x_axis,EE_EE_t_S1_global, EE_EE_t_CI_S1_global(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,EE_EE_t_S3_global, EE_EE_t_CI_S3_global(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,EE_EE_t_local, EE_EE_t_CI_local(:,2),'-.og','LineWidth',1),
    hold on,
    legend('EE JT-CoMP-NOMA','EE Conventional NOMA','ILO');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Average energy efficiency (b/s/Joule)');
    title(sprintf('EE - Optimized with %s (rho = %0.2f, kappa = %0.2f, P_{fix} = %0.2f) evaluated with PCM 1 (rho = %0.2f, kappa = %0.2f, P_{fix} = %0.2f)',PCM_text, data.rho, data.kappa, data.P_fix,0, 0, 0));
    hAx.YScale='log';
%     end
    
    % JT only
    [EE_EE_t_S1_global_PCM_prop, EE_EE_t_CI_S1_global_PCM_prop] = mean_confidence_interval(R_tot_EE_S1_global(:,feasible_JT)./Pi_sys_EE_S1_global_prop(:,feasible_JT));
    [EE_EE_t_S1_global_PCM2, EE_EE_t_CI_S1_global_PCM2] = mean_confidence_interval(R_tot_EE_S1_global(:,feasible_JT)./Pi_sys_EE_S1_global_PCM2(:,feasible_JT));
    [EE_EE_t_S1_global_PCM1, EE_EE_t_CI_S1_global_PCM1] = mean_confidence_interval(R_tot_EE_S1_global(:,feasible_JT)./Pi_sys_EE_S1_global_PCM1(:,feasible_JT));
    
    figure;
    hAx=axes;
    errorbar(x_axis,EE_EE_t_S1_global_PCM1, EE_EE_t_CI_S1_global_PCM1(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,EE_EE_t_S1_global_PCM2, EE_EE_t_CI_S1_global_PCM2(:,2),'-.ok','LineWidth',1)
    hold on,
    errorbar(x_axis,EE_EE_t_S1_global_PCM_prop, EE_EE_t_CI_S1_global_PCM_prop(:,2),'-.og','LineWidth',1)
    hold on,
    legend('PCM 1', 'PCM 2', 'PCM Proposed');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Average energy efficiency (b/s/Joule)');
    %title('EE JT-CoMP-NOMA');
    title(sprintf('EE JT-CoMP-NOMA - Optimized with %s (rho = %0.2f, kappa = %0.2f, P_{fix} = %0.2f)',PCM_text, data.rho, data.kappa, data.P_fix));
    hAx.YScale='log';
end


% --- Energy consumption plots -----
if(plot_EC)
    [EC_EE_t_S1_global, EC_EE_t_CI_S1_global] = mean_confidence_interval(Pi_sys_EE_S1_global_PCM1(:,feasible_all));
    [EC_EE_t_S3_global, EC_EE_t_CI_S3_global] = mean_confidence_interval(Pi_sys_EE_S3_global_PCM1(:,feasible_all));
    [EC_EE_t_local, EC_EE_t_CI_local] = mean_confidence_interval(Pi_sys_EE_local_PCM1(:,feasible_all));
    
    figure;
    hAx=axes;
    errorbar(x_axis,EC_EE_t_S1_global,EC_EE_t_CI_S1_global(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,EC_EE_t_S3_global,EC_EE_t_CI_S3_global(:,2),'-.ok','LineWidth',1),
    hold on,
    errorbar(x_axis,EC_EE_t_local,EC_EE_t_CI_local(:,2),'-.og','LineWidth',1),
    hold on,
    plot(x_axis,ones(size(x_axis))*data.Pt*data.N_BSs)
    legend('EE JT-CoMP-NOMA','EE Conventional NOMA','ILO');
    xlabel('Minimum data rate requirement (Kbps)'); 
    ylabel('Average system energy consumption'); % TODO: is it in Joule or Watt?
    title('Global Solution');
    
    [EC_EE_t_S1_global, EC_EE_t_CI_S1_global] = mean_confidence_interval(Pi_sys_EE_S1_global_PCM1(:, feasible_JT & ~feasible_NOMA ));
    figure;
    hAx=axes;
    errorbar(x_axis,EC_EE_t_S1_global,EC_EE_t_CI_S1_global(:,2),'-+b','LineWidth',1)
    hold on,
    plot(x_axis,ones(size(x_axis))*data.Pt*data.N_BSs)
    legend('EE JT-CoMP-NOMA','EE Conventional NOMA','ILO');
    xlabel('Minimum data rate requirement (Kbps)'); 
    ylabel('Average system energy consumption'); % TODO: is it in Joule or Watt?
    title('Global Solution');
end

plot_individual_TP = false;
if (plot_individual_TP)
% samples x N_users
    R_idx = 4;
    Ri_JT = R_EE_S1_global(R_idx,:,1, feasible_JT & feasible_NOMA);
    Ri_NOMA = R_EE_S3_global(R_idx,:,1, feasible_JT & feasible_NOMA);
    Ri_JT = [squeeze(Ri_JT).', squeeze(R_EE_S1_global(R_idx,1:data.N_inner_users,2, feasible_JT & feasible_NOMA)).'];
    Ri_NOMA = [squeeze(Ri_NOMA).', squeeze(R_EE_S3_global(R_idx,1:data.N_inner_users,2, feasible_JT & feasible_NOMA)).'];
    
    data = {Ri_JT, Ri_NOMA};
    figure,boxplotGroup(data, 'PrimaryLabels', {'JT' 'Conv.'}, ...
      'SecondaryLabels',cellstr(["User 1","User 2", "edge user", "User 1 (BS2)","User 2 (BS2)"]), 'InterGroupSpace', 1, ...
       'GroupType','withinGroups')
    title(sprintf('Throughput - %s',PCM_text));
    ylabel('User data rate (b/s)')
end

if (plot_individual_EE)
    % --- individual EE plots -----
    N_users = data.N_inner_users+data.N_JT_users;
    
    tot_user_power_i_bs_S1 = zeros(length(x_axis),N_users,data.N_BSs,N_samples);
    tot_user_power_i_bs_S3 = zeros(length(x_axis),N_users,data.N_BSs,N_samples);
    for samp = 1:N_samples
        for x_idx = 1:length(x_axis)
            %for ii = 1:N_users
                for bs=1:data.N_BSs
                    Pib_EE_S1_global = Pvec2mat(data.gamma, true, data.Pi_EE_S1_global(x_idx,:,samp));
                    Pib_EE_S3_global = Pvec2mat(data.gamma, false, data.Pi_EE_S3_global(x_idx,:,samp));
                    [~, tot_user_power_i_bs_S1(x_idx,:,bs,samp)] = user_power_consumption_2(Pib_EE_S1_global(:,bs), data.gamma, data.rho, data.kappa, false);
                    [~, tot_user_power_i_bs_S3(x_idx,:,bs,samp)] = user_power_consumption_2(Pib_EE_S3_global(:,bs), data.gamma, data.rho, data.kappa, false);
                end
            %end
        end
    end
    
    % For each cell-center user
    for bs=1:data.N_BSs
        for ii = 1:N_users-1
            [EE_Ri_EE_S1, EE_Ri_EE_CI_S1] = mean_confidence_interval( squeeze(R_EE_S1_global(:,ii,bs,feasible_all)./tot_user_power_i_bs_S1(:,ii,bs,feasible_all)) );
            [EE_Ri_EE_S3, EE_Ri_EE_CI_S3] = mean_confidence_interval( squeeze(R_EE_S3_global(:,ii,bs,feasible_all)./tot_user_power_i_bs_S3(:,ii,bs,feasible_all)) );
            
            
            figure;
            hAx=axes;
            errorbar(x_axis,EE_Ri_EE_S1,EE_Ri_EE_CI_S1(:,2),'-+b','LineWidth',1)
            hold on,
            errorbar(x_axis,EE_Ri_EE_S3,EE_Ri_EE_CI_S3(:,2),'-.ok','LineWidth',1),
            hold on,
            legend('EE JT-CoMP-NOMA','EE Conventional NOMA');
            xlabel('Minimum data rate requirement (Kbps)');
            ylabel('Average energy efficiency (b/s/Joule)');
            hAx.YScale='log';
            title(sprintf('User %d BS %d - %s',ii,bs,data.PCM));
        end
    end
    
    [EE_Ri_EE_S1, EE_Ri_EE_CI_S1] = mean_confidence_interval( squeeze(R_EE_S1_global(:,N_users,1,feasible_all)./sum(tot_user_power_i_bs_S1(:,N_users,:,feasible_all),3)) );
	[EE_Ri_EE_S3, EE_Ri_EE_CI_S3] = mean_confidence_interval( squeeze(R_EE_S3_global(:,N_users,1,feasible_all)./tot_user_power_i_bs_S3(:,N_users,1,feasible_all) ) );
    figure;
    hAx=axes;
    errorbar(x_axis,EE_Ri_EE_S1,EE_Ri_EE_CI_S1(:,2),'-+b','LineWidth',1)
    hold on,
    errorbar(x_axis,EE_Ri_EE_S3,EE_Ri_EE_CI_S3(:,2),'-.ok','LineWidth',1),
    hold on,
    legend('EE JT-CoMP-NOMA','EE Conventional NOMA');
    xlabel('Minimum data rate requirement (Kbps)');
    ylabel('Average energy efficiency (b/s/Joule)');
    hAx.YScale='log';
    title(sprintf('User %d (Edge User) - %s',N_users,data.PCM));
    
end