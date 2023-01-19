function [R_JT, R_NOMA, R_LOCAL, Pi_sys_JT, Pi_sys_NOMA, Pi_sys_local, feasible_JT, feasible_NOMA, feasible_LOCAL] = prepare_data_box(data, kappa, rho, P_fix)
R_EE_S1_global_00 = data.R_EE_S1_global;
R_EE_S3_global_00 = data.R_EE_S3_global;
R_EE_local_00 = data.R_EE_local;
R_EE_S1_global_00(isnan(R_EE_S1_global_00)) = 0;
R_EE_S3_global_00(isnan(R_EE_S3_global_00)) = 0;
R_EE_local_00(isnan(R_EE_local_00)) = 0;
R_JT = squeeze(sum(sum(R_EE_S1_global_00,2),3));
R_NOMA = squeeze(sum(sum(R_EE_S3_global_00,2),3));
R_LOCAL = squeeze(sum(sum(R_EE_local_00,2),3));


if(data.s~=data.N_samples)
    interval = 1:data.s-1;
else
    interval = 1:data.s;
end
simu_samps_00 = zeros(1,data.N_samples);
simu_samps_00(interval) = 1;

feasible_JT = (~sum((data.Exit_EE_S1_global < 0),1)>0) & simu_samps_00 & sum(isnan(data.Exit_EE_S1_global),1)==0;
feasible_NOMA = (~sum((data.Exit_EE_S3_global < 0),1)>0) & simu_samps_00 & sum(isnan(data.Exit_EE_S3_global),1)==0;
feasible_LOCAL = (~sum((data.Exit_EE_local < 0),1)>0) & simu_samps_00 & sum(isnan(data.Exit_EE_local),1)==0;


x_axis = data.x_axis;
N_samples = data.N_samples;

% Computes the system power expenditure

if(exist('rho', 'var'))
    rho_aux = rho;
else
    rho_aux = data.rho;
end
if(exist('kappa','var'))
    kappa_aux = kappa;
else
    kappa_aux = data.kappa;
end
if(exist('P_fix','var'))
    P_fix_aux = P_fix;
else
    P_fix_aux = data.P_fix;
end

Pi_sys_JT = zeros(length(x_axis),N_samples);
Pi_sys_NOMA = zeros(length(x_axis),N_samples);
Pi_sys_local = zeros(length(x_axis),N_samples);
for r=1:length(data.R_Kbps)
    for ss=1:N_samples
        Pi_sys_JT(r,ss) = system_power_consumption(data.Pi_EE_S1_global(r,:,ss), data.gamma, rho_aux, P_fix_aux, kappa_aux, true, false);
        Pi_sys_NOMA(r,ss) = system_power_consumption(data.Pi_EE_S3_global(r,:,ss), data.gamma, rho_aux, P_fix_aux, kappa_aux, false, false);
        Pi_sys_local(r,ss) = system_power_consumption(data.Pi_EE_local(r,:,ss), data.gamma, rho_aux, P_fix_aux, kappa_aux, false, false);
    end
end 

end