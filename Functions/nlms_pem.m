function [f_hat,d,a_hat] = nlms_pem(F_hat,y,u,d,nA,mu_f,alpha,R,varargin)
Uk = fft(u);
Y = fft([zeros(R/2,1);y],R);
y_hat = ifft(Uk.*F_hat,'symmetric');
d_temp = y - y_hat(R/2+1:end);
% AR estimation 
[rd, lg]= xcorr([d_temp;d],'biased'); % autocorrelation function for AR model
rd(lg<0)=[]; % taking only positive lags
%     a_hat = levinsondurbin(rd,nA);  
if length(varargin)==1
    a_hat = varargin{1};
elseif length(varargin)==2
    if varargin{1}==0
        a_hat = levinson(rd,nA);
    else
        a_hat = varargin{2};
    end
else
    a_hat = levinson(rd,nA);
end
d = d_temp;
A_hat = fft(a_hat.',R);
if iscolumn(A_hat)~=1 % if A_hat is not a column vector, make it so. 
    A_hat = A_hat.';
end
tempE = ifft(A_hat.*(Y - Uk.*F_hat),'symmetric');
Epsilon = fft([zeros(R/2+nA,1); tempE(end-R/2+nA+1:end,1)],R);
mu = 1./(alpha + sum(abs(Uk.*A_hat).^2,2));
DeltaF_hat = mu.*conj(A_hat.*Uk).*Epsilon;
temp_DF = ifft(DeltaF_hat,'symmetric');
% temp_DF = [zeros(638,1);temp_DF(1:R/2-638,1)];
F_hat = F_hat + mu_f*fft(temp_DF(1:R/2,:,:),R);
%     F_hat = F_hat + mu_f*fft((eye(M) - Q*Q')*ifft(DeltaF_hat));
f_hat = ifft(F_hat,'symmetric');
f_hat = f_hat(1:R/2,1);
% f_hat = [zeros(638,1);f_hat(1:R/2-638,1)];
end