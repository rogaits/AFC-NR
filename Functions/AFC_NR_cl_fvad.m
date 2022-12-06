function output = AFC_NR_cl_fvad(filename,iSNR,plotting_flag,alpha,mu_f,nA,MetricsToCompute,R,asg_flag,lambda)
% This functions performs multi-channel acoustic feedback cancellation followed by
% noise reduction for a multi-microphone single-loudspeaker (MMSL) scenario. 
% It uses an M-channel PEM-based AFC algorithm followed by a rank-1
% M-channel NR stage. For more details refer to the following paper. 
% 
% Ruiz, S., van Waterschoot, T., Moonen M. Cascade algorithms for combined
% acoustic feedback cancelaltion and noise reduction. 2022. 
% 
% INPUT PARAMETERS
%   filename        = string with filename where the scenarios is saved
%   iSNR            = input SNR
%   plotting_flag   = 1 for showing plots 0 otherwise
%   alpha           = Regularization parameter for the PEM-based AFC algorithm
%   mu_f            = Step size for the PEM-based AFC algorithm
%   nA              = Order of the AR model for the PEM-based AFC algorithm
%   MetricsToCompute= Cell array with string inside the cells specifying
%   the metrics to compute e.g. {'stoi','sd'}
%   R               = Window size
%   asg_flag        = 1 for computing ASG with MWF filters, 0 otherwise
%   lambda          = Forgetting factor for covariance estimation
% 
% OUTPUT PARAMETERS
%   output= Cell array containing  the following variables:
%       t2                  = time vector
%       ASG_dB_pem          = ASG values for each time frame
%       Mis                 = Mis values for each time frame
%       new_metrics_table   = Computed metrics on the estimated signal
%       nmic_metrics_table  = Computed metrics on the microphone signal
%       estimation          = Desired signal estimate
%       Ggain_profile       = Forward path gain profile

% Author: Santiago Ruiz
% Date: November 2022
% STADIUS - ESAT 
% KU Leuven, Belgium 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Function starts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('../Scenarios'))
addpath(genpath('../Metrics'))
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           DATA GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(filename)
idx_ref = 1;
delay_nr = R;       % Delay for the NR stage in samples
delay_afc = R/2;    % Delay for the AFC stage in samples
Gdelay = 10;        % Extra delay to account for other processing stages
delay_fwd = Gdelay+delay_nr+delay_afc; % Forward path delay NR+AFC stages
%% Initial ASG computation 
if asg_flag==0
    if R ~= nF
       for i=1:M
        [msg_dB(i),omega_msg(i),idx_max_msg(i)] = msg_calculation([zeros(delay_fwd,1);fTrue(:,i)],dfilt.df2sos([1 0 0 1 0 0]),0,nF);
        end
    end
else
    H_mwf = [ones(R/2+1,1) zeros(R/2+1,M-1)];
    [msg_dB,omega_msg,idx_max_msg] = msg_calculation([zeros(delay_fwd,M);fTrue],dfilt.df2sos([1 0 0 1 0 0]),0,R/2+1,H_mwf);
end
%% Gain profile
K_msg = min(msg_dB);
L1 = length(Lf:R/2:Nt);
K1 = 10.^((K_msg-15)/20);
K2 = 10.^((K_msg)/20);
Ggain_profile = [linspace(K1,K1,ceil(L1/2)) linspace(K1,K2,ceil(L1/4)) linspace(K2,K2,ceil(L1/4))];
Ggain_profile = Ggain_profile(1:L1);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Initialization of variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimation = zeros(Nt+Lf+Gdelay,1);     % desired signal estimation
win = sqrt(hann(R,'periodic'));         % WOLA analysis and synthesis window 
% Signal covariance mattrices initializations
dimR = M;                                % Covariance signal dimension
Ryy = zeros(dimR,dimR,R/2+1);            % Speech-plus-noise covariance matrix
Rnn = repmat(1e-6*eye(dimR),1,1,R/2+1);  % noise-only covariance matrix 
if ~exist('lambda')
    lambda = 0.9997;                        % forgetting factor
end
Ryy_initialized = zeros(R,1);       % FLAG for Ryy initialization per frequency bin
Rnn_initialized = zeros(R,1);       % FLAG for Rnn initialization per frequency bin
Ryy_samples = zeros(R,1);           % Number of samples used for estimating Ryy 
Rnn_samples = zeros(R,1);           % Number of samples used for estimating Rnn
W = [eye(1); zeros(M-1,1)];         % Initial MWF filter
l = 1;                              % frame index
% y_hat_nr = zeros(Nt+Lf+Gdelay+delay_nr,1);     % desired signal estimation after NR stage
% u_hat_nr = zeros(Nt+Lf+Gdelay+delay_nr,1);     % desired signal estimation after NR stage
d_hat = zeros(Nt+Lf+Gdelay+delay_nr,1);     % desired signal estimation

F_hat = zeros(R,M);     % Estimated feedback paths in the frequency domain
f_hat = zeros(R/2,M);   % Estimated feedback paths in the time domain
d = zeros(R/2,M);       % Desired signals blocks in the afc stage
% Frequency - domain VAD
alpha_vad = db2pow(0); % 3dB is the default value in the function
myvad_freq = vad_freq(v0, Fs, win, R, R/2,alpha_vad);

fTrue(R,1) = 0; 
if Lf>=R   
    k_init = Lf;
else
    k_init = R;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Closed-loop data generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = k_init:R/2:Nt
    Ggain = Ggain_profile(l);       % Forward path gain based in Gain profile
    idx_f2 = k-Lf+1:k;              % indexes for signal generation
    xf2(idx_f2,:) = fftfilt(fTrue,u2(idx_f2,:));                % Feedback component
    x2(idx_f2,:) = v(idx_f2,:) + xf2(idx_f2,:) + n(idx_f2,:) ;  % microphone signals
    y2(idx_f2,:) = [u2(idx_f2,:) x2(idx_f2,:)];                 % Nx1 Signal vector 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AFC stage
    idx_f = k:-1:k-R+1;             % OLS indexes
    % indexes for microphone and loudspeaker signals
    idx_u = fliplr(idx_f);          % Loudspeaker indexes
    idx_y = idx_u(end-R/2+1:end);   % Microphone indexes
    % Filter estimation using PEM-based AFC algorithm 
    if vad_speaker(k) == 1 
        % testing: using known AR model coefficients
        if var(y2(idx_u,1))>0
            for m=1:M
                [f_hat(:,m),d(:,m),a_hat(m,:)] = nlms_pem(F_hat(:,m),y2(idx_y,L+m),y2(idx_u,1),d(:,m),nA,mu_f,alpha,R);
            end
        else
            a_hat_aux = [1 zeros(1,nA)]; % this is done to avoid NaN values
%           in the AR estimation because the loudspeaker signal is zero
            for m=1:M
                [f_hat(:,m),d(:,m),a_hat(m,:)] = nlms_pem(F_hat(:,m),y2(idx_y,L+m),y2(idx_u,1),d(:,m),nA,mu_f,alpha,R,a_hat_aux);
            end
        end
        x_hat(idx_y+delay_afc,:) = d; %storing the estimated signal
    else
        % When there is no signal in the loudspeaker, signals are filtered
        % with previous estimates
        U2 = fft(y2(idx_u,1));
        y_hat = ifft(U2.*F_hat,'symmetric');
        x_hat(idx_y+delay_afc,:) = y2(idx_y,L+1:N) - y_hat(R/2+1:end,:);    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Computation of some AFC metrics
    mse(l,:) = mean((fTrue(1:R/2,:) - f_hat).^2);
    Mis(l,:) = 20*log10(vecnorm(fTrue(1:R/2,:) - f_hat,2,1)./vecnorm(fTrue(1:R/2,:),2,1));
    F_hat = fft(f_hat,R);
    F_true = fft(fTrue,R);
    if asg_flag==0
        for m=1:M
            ASG_dB_pem(l,m) = asg_computation(fTrue(:,m),f_hat(:,m),Gdelay,R,msg_dB(m),idx_max_msg(m));
        end
    else        
        ASG_dB_pem(l) = asg_computation(fTrue,f_hat,delay_fwd,R,msg_dB,idx_max_msg,H_mwf);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idx_wola = k-R+1:k;                 % indexes for WOLA
    Y = fft(win.*x_hat(idx_wola,:)).';  % Microphone signal frame
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NR stage
    for kbin=1:R/2+1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % VAD check of speech activity
        if myvad_freq(kbin,l) == 1   % speech is active
            Ryy_samples(kbin) = Ryy_samples(kbin)+1;
            if Ryy_samples >= N+1
                Ryy_initialized(kbin) = 1;
            end
        else
            Rnn_samples(kbin) = Rnn_samples(kbin)+1;
            if Rnn_samples >= N+1
                Rnn_initialized(kbin) = 1;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Ryy and Rnn estimation based on VAD
        Y_aux(:,kbin) = Y(:,kbin);
        if myvad_freq(kbin,l) == 1
            Ryy(:,:,kbin) = lambda*Ryy(:,:,kbin) + (1-lambda)*(Y_aux(:,kbin)*Y_aux(:,kbin)');
        else
            Rnn(:,:,kbin) = lambda*Rnn(:,:,kbin) + (1-lambda)*(Y_aux(:,kbin)*Y_aux(:,kbin)');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GEVD-based diagonalization to estimate MWF 
        if Ryy_initialized(kbin) == 1 && Rnn_initialized == 1 
            W = mwf_GEVD2(Ryy(:,:,kbin),Rnn(:,:,kbin),1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Y_hat(kbin,:) = W(:,1)'*Y_aux(:,kbin);  % MWF filtering for microphone signals only
        if asg_flag==1
            H_aux= W';
            % selecting the filter coefficients that correspond to the mics
            H_mwf(kbin,:) = H_aux(1,1:M);   
        end
    end
    % WOLA synthesis
    [d_hat,d_u] = stft_synthesis(l,d_hat,R,Y_hat,0,win); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % looudspeaker signal generation
    idx_futuresamples2 = idx_wola+delay_fwd;            % Current samples are delayed 
    u2(idx_futuresamples2,:) = Ggain*d_hat(idx_wola,1); % forward paht gain applied
    l=l+1;          % Frame index increase
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     END of Closed-loop data generation
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Metrics computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
estimation = d_hat(delay_afc+1:end-2*R); % compensating for the delay for the metrics
v = v(1:length(estimation),:);           % Speech component in microphone signals
x2 = x2(1:length(estimation),1);         % Microphone signals
vad_speaker = vad_speaker(1:length(estimation));    % Time-domain VAD used for metrics computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       SD, STOI, iSNR, SFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input SNR check
iSNR_ = 10*log10(var(v(vad_speaker==1,:))./var(n))
% Signal-to-feedback ratio computation
SFR = 10*log10(var(v(vad_speaker==1,:))./var(xf2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Metrics table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[new_metrics_table, varNames] = Metrics_Compute(v(vad_speaker==1,1), estimation(vad_speaker==1,1),Fs,MetricsToCompute);
[nmic_metrics_table, varNames] = Metrics_Compute(v(vad_speaker==1,1), x2(vad_speaker==1,1),Fs,MetricsToCompute);
new_metrics_table
nmic_metrics_table
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotting_flag==1
    M2 = 512; % window size for plotting
    figure('Name',['AFC-NR: Spectrograms iSNR=' num2str(iSNR)])
    subplot(411)
    spectrogram(v(:,1),hann(M2),M2/2,M2,Fs,'yaxis')
    title('Near-end')
    caxis([-150 0])
    subplot(412)
    spectrogram(x2(:,1),hann(M2),M2/2,M2,Fs,'yaxis')
    title('Microphone signal')
    caxis([-150 0])
    subplot(413)
    spectrogram(u2(:,1),hann(M2),M2/2,M2,Fs,'yaxis')
    title('Loudspeaker signal')
    caxis([-150 0])
    subplot(414)
    spectrogram(estimation(:,1),hann(M2),M2/2,M2,Fs,'yaxis')
    title('Estimated signal')
    caxis([-150 0])
    annotation('textbox',[.9 .3 .4 .5],'String','my text','FitBoxToText','on')
    annotation('textbox', [0, 0.5, 0, 0], 'string', 'My Text')

    figure('Name',['AFC-NR: Time-domain signals iSNR=' num2str(iSNR)]);hold on;grid on
    t = (0:length(x2(:,1))-1)/Fs;
    plot(t,x2(:,1),'DisplayName','Microphone signal')
    plot(t,v(:,1),'DisplayName','Desired speech component at reference microphone')
    plot(t,estimation(:,1),'DisplayName','Estimated signal')
    plot(t,vad_speaker)
    legend('Fontsize',14)
    xlabel('Time [s]','Fontsize',14)
    ylabel('Amplitude','Fontsize',14)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Mis, ASG, MSE plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t2 = (Lf:R/2:Nt)/Fs;
    figure('Name',['AFC-NR: Mis, ASG, MSE iSNR=' num2str(iSNR)])
    subplot(311)
    plot(t2,ASG_dB_pem,'k','LineWidth',2)
    grid on;ylim([0 50]);grid minor;xlim([0 max(t2)]);
    title('ASG')
    xlabel('Time (s)')
    subplot(312)
    plot(t2,Mis,'k','LineWidth',2)
    grid on;grid minor;xlim([0 max(t2)]);
    title('Misalignment')
    xlabel('Time (s)')
    subplot(313)
    semilogy(t2,mse,'k','LineWidth',2)
    grid on;grid minor;xlim([0 max(t2)]);
    title('MSE')
    xlabel('Time (s)')
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   \hat{f} and f_true plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    idx_ref = 1;
    figure('Name',['AFC-NR: Estimated and True feedback path iSNR=' num2str(iSNR)])
    plot(f_hat(:,idx_ref),'Displayname','$\hat{f}_r$','LineWidth',1)
    hold on;grid on;xlim([0 length(fTrue(:,idx_ref))])
    plot(fTrue(:,idx_ref),'--','Displayname','$f_r$','LineWidth',1)
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    legend('Interpreter','Latex','FontSize',20)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output = {t2,ASG_dB_pem,Mis,new_metrics_table,nmic_metrics_table,estimation,Ggain_profile};
end
function [estimation,d] = stft_synthesis(l,estimation,R,Y,delay,win)
% l = current frame
% estimation = concatenated signal
% R window size
% Y signal frame to synthesize
% STFT synthesis Batch
    if ~iscolumn(Y)
        Yest = Y.';
    else
        Yest = Y;
    end
    blockest = real(ifft([Yest;conj(flipud(Yest(2:R/2,:)))],'symmetric'));
    idxs = ((l-1)*R/2+1:(l+1)*R/2) + delay;
    estimation(idxs,:)=estimation(idxs,:)+win.*blockest;
    d = sqrt(hann(R)).*blockest;
end
