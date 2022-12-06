function SD = SD_measure_fbcancellation(v,d,fs,Nsec)
% d = processed signal
% v = clean signal
% Nsec = signal duration in seconds

% % switch fs
% %     case 8000
% %         P_psd = 1024;
% %     case 16000
% %         P_psd = 2048;
% %     case 44100
% %         P_psd = 4096;
% % end

P_psd = 512;

[~,~,~,Pd] = spectrogram(d,P_psd,P_psd/2,P_psd,fs,'yaxis');

[~,F,T,Pv] = spectrogram(v,P_psd,P_psd/2,P_psd,fs,'yaxis');

w_ERB_Ann = Ann_ERBweightGISO(F);

LSD_indices = find(T-Nsec/2-(P_psd/(2*fs)) >= 0);
SD = zeros(length(LSD_indices),1);
for k = 1:length(LSD_indices),
    LSD_index = LSD_indices(k);
    squared_log_spectral_ratio = (10*log10(Pd(:,LSD_index)./Pv(:,LSD_index))).^2;
    SD(k) = sqrt(w_ERB_Ann'*squared_log_spectral_ratio);
end

% Mean and max SD
SD = [mean(SD) max(SD)];

end
