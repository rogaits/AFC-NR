clearvars;clc;close all;
addpath(genpath('Scenarios'))
addpath(genpath('Functions'))
addpath(genpath('Metrics'))

plotting_flag = 1;
iSNR = 0;
filename = ['Scenarios/data_cl_cas_ARmodels_Hm_SNR' num2str(iSNR)];
metrics2compute = {'stoi','sd'};
nA = 30;
R = 1024; % length of IRs is 1024 samples
if iSNR==10
    alpha = 90e0;
    mu_f = 0.4;
elseif iSNR==20
    alpha = 5e0;
    mu_f =3e-1;
else
    alpha = 90e0;
    mu_f = 0.4;
end
asg_flag = 1; % 1 -> computes ASG considering the MWF filters
out_rank2_NR_AFC = rank2_NR_AFC_cl_fvad(filename,iSNR,plotting_flag,alpha,mu_f,nA,metrics2compute,R,asg_flag);
out_rank1_NR_AFC = rank1_NR_AFC_cl_fvad(filename,iSNR,plotting_flag,alpha,mu_f,nA,metrics2compute,R,asg_flag);
out_AFC_NR = AFC_NR_cl_fvad(filename,iSNR,plotting_flag,alpha,mu_f,nA,metrics2compute,R,asg_flag);
% scenario_params = load(filename)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting ASG, Mis for all algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=out_rank2_NR_AFC{1};
figure('Name',['ASG and Mis, iSNR=' num2str(iSNR)])
subplot(211)
plot(t,out_rank2_NR_AFC{2},'k','LineWidth',2)
grid on;hold on;ylim([-5 40]);grid minor;xlim([0 max(t)]);
plot(t,out_rank1_NR_AFC{2},'r--','LineWidth',2)
plot(t,out_AFC_NR{2},'b','LineWidth',2)
legend({'Rank-2 NR-AFC','Rank-1 NR-AFC','AFC-NR',})
title('ASG')
xlabel('Time (s)')
subplot(212)
plot(t,out_rank2_NR_AFC{3},'k','LineWidth',2)
grid on;hold on;grid minor;xlim([0 max(t)]);
plot(t,out_rank1_NR_AFC{3},'r--','LineWidth',2)
plot(t,out_AFC_NR{3}(:,1),'b','LineWidth',2)
% plot(t,out_NR_AFC1_2{3},'k','LineWidth',2)
title('Misalignment')
xlabel('Time (s)')
set(findall(gcf,'-property','FontSize'),'FontSize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   metrics table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metrics = [out_rank2_NR_AFC{4};out_rank1_NR_AFC{4};out_AFC_NR{4}];
algo_name = {'Rank-2 NR-AFC';'Rank-1 NR-AFC';'AFC-NR'};
metrics = [table(algo_name,'VariableNames',{'Algorithm'}) metrics]
