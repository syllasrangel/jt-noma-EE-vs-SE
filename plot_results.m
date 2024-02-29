function [] = plot_results(filepath, plot_outage, plot_iter, plot_throughput, plot_EE, plot_EC, plot_individual_EE)

data = load(filepath);

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

feasible_and_same_transmiting_BS = feasible_JT & feasible_NOMA & (squeeze(data.Pi_EE_S1_global(1,end,:)) == 0).';
%feasible_all = feasible_and_same_transmiting_BS;


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
    %feasible_and_same_transmiting_BS = feasible_JT & feasible_NOMA & (squeeze(data.Pi_EE_S1_global(1,6,:)) == 0).';
    %feasible_all = feasible_and_same_transmiting_BS;
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

plot_individual_TP = true;
if (plot_individual_TP)
% samples x N_users
    R_idx = 4;
    Ri_JT = R_EE_S1_global(R_idx,:,1, feasible_JT & feasible_NOMA);
    Ri_NOMA = R_EE_S3_global(R_idx,:,1, feasible_JT & feasible_NOMA);
    Ri_JT = [squeeze(Ri_JT).', squeeze(R_EE_S1_global(R_idx,1:data.N_inner_users,2, feasible_JT & feasible_NOMA)).'];
    Ri_NOMA = [squeeze(Ri_NOMA).', squeeze(R_EE_S3_global(R_idx,1:data.N_inner_users,2, feasible_JT & feasible_NOMA)).'];
    
    box_data = {Ri_JT, Ri_NOMA};
    figure,boxplotGroup(box_data, 'PrimaryLabels', {'JT' 'Conv.'}, ...
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