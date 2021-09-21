clear
% close all
clc

addpath('F:\Research\MATLAB codes')
addpath('F:\Research\MATLAB')
addpath('F:\Research\MATLAB\Spectrogram\')

wav_files = dir('*_ramp.wav');
wav_files = {wav_files.name};

pairs = ['1 vs. 4'; '2 vs. 5'; '3 vs. 6'; '4 vs. 7'; '5 vs. 8'; '6 vs. 9'];

for i = 1 : length(wav_files)

    file_name = char( wav_files(i) );
    if(contains(file_name, 'vs')), continue; end
    
    if(contains(file_name, 'chirp'))
        title_ = 'chirp';
    else
        title_ = 'speech';
    end
    
    indx = str2double(file_name(strfind(file_name, '_')+1:end-length('_ramp.wav')));
    
    disp(file_name)
    [d, fs] = audioread( file_name );
    R = rms(d);
    if(R(1)>R(2))        
        d = d(:,1);
    else
        d = d(:,2);
    end
    sp(d, fs, 1024, 160, 80);
    
    set(gca, 'FontSize', 15);
    title( [title_ '__slopes ' pairs(indx,:) ]  )
    
    temp = [title_ '__slopes ' pairs(indx,:) ];
    temp = strrep( temp, '.', '' );    
    
    saveas(gcf, [temp '_spec.jpg'])
    audiowrite( [temp '_ramp_5ms.wav'], d, fs );
    
    plot( 1000*(1:length(d))/fs, d )   
    xlabel('Time (ms)')
    set(gca, 'FontSize', 15);
    title( [title_ '__slopes ' pairs(indx,:) ]  )
    axis tight
    saveas(gcf, [temp '_time.jpg'])
    
%     input('')
    
end
