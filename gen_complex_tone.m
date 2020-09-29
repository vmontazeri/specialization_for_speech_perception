function [cmplx_tone, trans_tone, stable_tone] = gen_complex_tone(f_start, f_end, f_stable, trans_dur, stable_dur, FS)
% *Function to generate a complex tone that consists of a stable part with 
% constant frequency as well as a tone with varying frequency. 
% f_start: Initial frequency of the varying frequency tone.
% f_end: Final frequency of the varying frequency tone.
% f_stable: Freuquency of the stable tone. 
% trans_dur: Duration of the varying-frequency tone in ms. 
% table_dur: Duration of the stable tone in ms. 
% FS: Sampling frequency. 
% Vahid Montazeri, 3/22/2020*

stable_tone = [];
trans_tone = [];

if(f_stable > 0)
    stable_tone = gen_tone(f_stable, stable_dur, FS, 1);
end
if(f_start * f_end > 0)
    trans_tone = freq_varying_tone(f_start, f_end, trans_dur, FS, 1);
end
cmplx_tone = [trans_tone(:); stable_tone(:)];
