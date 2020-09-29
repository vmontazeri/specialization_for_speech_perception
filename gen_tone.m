function tone = gen_tone(freq, L, FS, rnd)
% *Function to generate a tone
% freq: Desired frequency. 
% L: Tone duration in milliseconds. 
% FS: Sampling frequency. 
% rnd: if set 1, the tone will have random phase. Otherwise, phase will be
% zero. 
% Vahid Montazeri, 3/22/2020*

TWO_PI = 2 * pi;
time = TWO_PI * (1/FS:1/FS:L/1000);
if(rnd)
    tone = sin( time * freq + rand(1,1) );
else
    tone = sin( time * freq );
end