function vad_output = vad(x,th,win)

energy_threshold = th;
% win = w;

abs_x = abs(x);
signal_length = length(x);
vad3 = zeros(size(x));
for k=1 : signal_length
    if abs_x(k) > energy_threshold
        vad2(k) = 1;
    else 
        vad2(k) = 0;
    end
end
    
for k=win+1 : signal_length-win
    if sum(vad2(k-win:k+win)) > win/2
        vad3(k) = 1;
    else 
        vad3(k) = 0;
    end
end

vad_output = vad3;