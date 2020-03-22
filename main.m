clear
close all
clc

FS = 16e3;

freq = 750; 
L = 250;
tone1 = gen_tone(freq, L, FS, 1);
v_tone1 = freq_varying_tone(250, 750, 50, FS, 1);
tone1 = [v_tone1(:); tone1(:)];

freq = 1250; 
L = 250;
tone2 = gen_tone(freq, L, FS, 1);
v_tone2 = freq_varying_tone(1750, 1250, 50, FS, 1);
tone2 = [v_tone2(:); tone2(:)];

freq = 2500; 
L = 250;
tone3 = gen_tone(freq, L, FS, 1);
consonant = input('enter consonant ''d'' or ''g''\n', 's');
if( strcmpi(consonant, 'd') )
    disp('generating d')
    v_tone3 = freq_varying_tone(3250, 2500, 50, FS, 1);
else
    disp('generating g')
    v_tone3 = freq_varying_tone(1750, 2500, 50, FS, 1);
end
tone3 = [v_tone3(:); tone3(:)];

tone = ( tone1 + tone2 + tone3 ) / 3;


soundsc( [zeros(FS/2, 1); tone], FS )
spectrogram(tone, 128, 32, 1024, 'yaxis')
% imagesc(abs(fft(tone2')))