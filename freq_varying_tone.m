function tone = freq_varying_tone(freq0, freq1, L, FS, rnd)
% *Function to generate a tone with varying frequency 
% freq0: Initial frequency 
% freq1: Final frequency.  
% FS: Sampling frequecy. 
% rnd: if set 1, the tone will have random phase. Otherwise, phase will be
% zero. 
% Vahid Montazeri, 3/22/2020*

TWO_PI = 2 * pi;
time = 1/FS:1/FS:L/1000;
f_in = linspace(freq0, freq1, length(time));
phase_in = cumsum(f_in/FS);
if(rnd)
    tone = sin( TWO_PI * phase_in + rand(1,1) );
else
    tone = sin( TWO_PI * phase_in );
end
