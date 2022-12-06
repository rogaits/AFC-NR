function [msg,omega_msg,index_max] = msg_calculation(loop_tf_num,loop_filterobject,phase_shift,n_freqz,varargin)
%msg_calculation Calculates the MSG of a closed-loop system based with a
%                FIR as feedback path and a Gain + delay as forward path
%
% 
%
% function [msg,omega_msg] = msg_calculation(loop_tf_num,...
%                            loop_filterobject,phase_shift,n_freqz)
%
% Inputs:
%    loop_tf_num        - open loop TF, e.g. '[zeros(delay,1);f]'
%    loop_filterobject  - pass the argument 'dfilt.df2sos([1 0 0 1 0 0])'
%    phase_shift        - phase shift
%    n_freqz            - number of coefficients used for the FFT
%
% Outputs:
%    msg          - maximum stable gain
%    omega_msg    - frequencies exhibiting possible instability
%
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% MATLAB version: R2017a
%
% See also: 
%
% Authors: Toon van Waterschoot, Giuliano Bernardi
% KU Leuven, Department of Electrical Engineering (ESAT/STADIUS)
% email: giuliano.bernardi@esat.kuleuven.be
% Website: http://homes.esat.kuleuven.be/~gbernard/index.html
% Created: 10-July-2009; Last revision: 12-November-2017
%
% This code is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your %
% option) any later version.
%
% This code is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see http://www.gnu.org/licenses/.

%------------- BEGIN CODE --------------

if isempty(varargin) ~= 1
    % This section computes the MSG when a beamforming operation is
    % performed in the closed-loop system. 
    H = varargin{1}; % STFT filter coefficients for MSG computation (MWF)
   
    M = size(loop_tf_num,2); % number of microphones
    for m=1:M
        [fTrue_freqz(:,m),omega] = freqz(loop_tf_num(:,m),1,n_freqz);
    end
    Y_H_fTrue = H.*fTrue_freqz; % multiply each microphone with its corresponding STFT filter

    loop_tf_freqz = sum(Y_H_fTrue,2); % sum over each microphone to reduce to a single channel MSG problem
    n = n_freqz; % desired order of the numerator
    m = 0; % desired order of the denominator
    [loop_tf_num,loop_tf_den] = invfreqz(loop_tf_freqz,omega,n,m);
    if length(varargin)>1
        F_hat = varargin{2};
        FG_gain_tf = loop_tf_freqz - F_hat;
    end
else
    [loop_tf_freqz,omega] = freqz(loop_tf_num,1,n_freqz);
    loop_tf_den = 1;
end
% [loop_filterobject_freqz,omega]    = freqz(loop_filterobject,n_freqz);
loop_freqz                         = loop_tf_freqz;%.*loop_filterobject_freqz;
loop_phase                         = (unwrap(angle(loop_freqz))+phase_shift)/(2*pi);
index_crossover_lower              = find(diff(ceil(loop_phase)));
omega_crossover_lower              = omega(index_crossover_lower);
omega_crossover_upper              = omega(index_crossover_lower+1);
loop_phase_crossover_lower         = loop_phase(index_crossover_lower);
loop_phase_crossover_upper         = loop_phase(index_crossover_lower+1);
omega_crossover_interpolated       = omega_crossover_upper -...
                                   abs(loop_phase_crossover_upper ...
                                   - ceil(loop_phase_crossover_upper)).*...
                                   ((omega_crossover_upper-omega_crossover_lower)./...
                                   abs(loop_phase_crossover_upper-loop_phase_crossover_lower));


                               
loop_tf_freqz_crossover            = freqz(loop_tf_num,loop_tf_den,omega_crossover_interpolated);
loop_filterobject_freqz_crossover  = freqz(loop_filterobject,omega_crossover_interpolated);
loop_freqz_crossover               = loop_tf_freqz_crossover.*loop_filterobject_freqz_crossover;

% figure;
% plot(omega,20*log10(abs(loop_freqz)));
% hold on;
% plot(omega_crossover_interpolated,20*log10(abs(loop_freqz_crossover)),'ro');
% figure;
% plot(omega,angle(loop_freqz));
% hold on;
% plot(omega_crossover_interpolated,angle(loop_freqz_crossover),'ro');
% figure;
% plot(omega,unwrap(angle(loop_freqz))/(2*pi));
% hold on;
% plot(omega_crossover_lower,loop_phase_crossover_lower,'ro');
% plot(omega_crossover_upper,loop_phase_crossover_upper,'g*');
% plot(omega_crossover_interpolated,loop_phase_crossover_upper,'bs');

[max_loop_gain_crossover,index_max] = max(abs(loop_freqz_crossover));
msg = -20*log10(max_loop_gain_crossover);
omega_msg = omega_crossover_interpolated(index_max);
