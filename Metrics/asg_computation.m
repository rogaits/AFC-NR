function ASG_dB = asg_computation(fTrue,f_hat,Gdelay,R,Kmsg_dB,idx_max,varargin)
if isempty(varargin) == 1
    if size(fTrue,1)>size(f_hat,1) 
        fr = fTrue(1:length(f_hat),1) - f_hat;
    elseif size(fTrue,1)==size(f_hat,1) 
        fr = fTrue - f_hat;
    else
        fTrue(R,1) = 0;
        fr = fTrue - f_hat;
    end
    [msg_dB,omega_msg,idx_max] = msg_calculation([zeros(Gdelay,1);fr],dfilt.df2sos([1 0 0 1 0 0]),0,R);
    % Ftrue = fft(fTrue,R);
    % ASG = abs(Ftrue(idx_max) - Fhat(idx_max));
    ASG_dB = msg_dB - Kmsg_dB;
elseif length(varargin) ==1 && isvector(f_hat)
    H = varargin{1};
    M = size(fTrue,2);
    if size(fTrue,1)>size(f_hat,1) 
        f_hat(size(fTrue,1),1)=0;
        fr = fTrue(1:length(f_hat),1) - f_hat;
    elseif size(fTrue,1)==size(f_hat,1) 
        fr = fTrue - f_hat;
    else
        fTrue(R,1) = 0;
        fr = fTrue - f_hat;
    end
%     f_open_loop = [zeros(Gdelay,M);[fr fTrue(1:R/2,2:end)]];
    f_open_loop = [zeros(Gdelay,M);[fr fTrue(:,2:end)]];
    F_hat = fft(f_hat,R);
    [msg_dB,omega_msg,idx_max] = msg_calculation(f_open_loop,dfilt.df2sos([1 0 0 1 0 0]),0,R/2+1,H,F_hat(1:R/2+1));
    % Ftrue = fft(fTrue,R);
    % ASG = abs(Ftrue(idx_max) - Fhat(idx_max));
    ASG_dB = msg_dB - Kmsg_dB;
elseif length(varargin) == 1 && ~isvector(f_hat)
    H = varargin{1};
    M = size(fTrue,2);
    if size(fTrue,1)>size(f_hat,1) 
        f_hat(size(fTrue,1),1)=0;
        fr = fTrue(1:length(f_hat),:) - f_hat; 
    elseif size(fTrue,1)==size(f_hat,1) 
        fr = fTrue - f_hat;
    else
        fTrue(R,1) = 0;
        fr = fTrue - f_hat;
    end
    f_open_loop = [zeros(Gdelay,M);fr];
    F_hat = fft(f_hat,R);
    [msg_dB,omega_msg,idx_max] = msg_calculation(f_open_loop,dfilt.df2sos([1 0 0 1 0 0]),0,R/2+1,H,F_hat(1:R/2+1,:));
    % Ftrue = fft(fTrue,R);
    % ASG = abs(Ftrue(idx_max) - Fhat(idx_max));
    ASG_dB = msg_dB - Kmsg_dB;
else
    disp('Error computing the ASG, wrong number of input parameters!')
end
% ftrue =  time-domain true feedback path
% Fhat = DFT of the estimated feedback path
% f = frequency bin indexes where instabilities occur
end