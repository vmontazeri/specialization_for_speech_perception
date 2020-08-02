% *code for experiments described in 'A specialization for speech perception'.
% Vahid Montazeri, 3/22/2010 '*
if(exist('conditions'))
    recovery = input('Recovery mode? (Y/N)\n', 's');
    if(strcmpi(recovery, 'Y'))
        recovery_mode = 1;
    else
        recovery_mode = 0;
    end
else
    recovery_mode = 0;
end

if(~recovery_mode)
    clear
    close all
    clc
    
    recovery_mode = 0;
    
    valid_answer = 0;
    while(~valid_answer)
        clc; disp('Enter mode (1) or (2):'); disp('(1) debug mode'); disp('(2) test mode');
        s = input('Enter your answer: ', 's');
        if(strcmpi(s, '1')), debug = 1; valid_answer = 1; elseif(strcmpi(s, '2')), debug = 0; valid_answer = 1; else, disp('invalid answer'); pause(1); end
    end
    
    listener_code = input('Enter listener code or press enter to generate an automatic code:\n', 's');
    listener_code = [listener_code '_' strrep(strrep(strrep(char(datetime), ':', '_'), '-', '_'), ' ', '_')];
end

% training( debug, listener_code );
F1_f_stable = 750;
F1_f_start    = 250;
F1_f_end     = 750;

F2_f_stable = 1250;
F2_f_start    = 1750;
F2_f_end     = 1250;

F3_f_stable          = 2500;
F3_min_trans_f   = 1750;
F3_max_trans_f  = 3250;

trans_dur = 50;
stable_dur = 250;
FS = 16e3;
left_ear = 1;
right_ear = 2;

silence = zeros(0.5*FS,2);
tone_len = (trans_dur + stable_dur)*FS/1000;
w = hanning(tone_len + 2*length(silence));

transition_indx = 1:trans_dur*FS/1000;
stable_indx = trans_dur*FS/1000+1: tone_len;

trial_per_level = 10;

%% experiment 1
factor1_levels = [1; 2; 3; 4];
factor2_levels = 1:9;
factor1_levels_1 = factor1_levels(randperm(length(factor1_levels)));

if(~recovery_mode)
    conditions = [];
    for i = 1 : length(factor1_levels)
        if(factor1_levels(i)==2)
            factor1_levels_2 = repelem(factor1_levels_1(i), trial_per_level);
            conditions = [conditions; [factor1_levels_2(:) -1*ones(trial_per_level,1)]];
        else
            factor2_levels_1 = factor2_levels(randperm(length(factor2_levels)));
            factor2_levels_2 = repelem(factor2_levels_1, trial_per_level);
            factor2_levels_3 = factor2_levels_2( randperm(length(factor2_levels_2)) );
            factor1_levels_2 = repelem(factor1_levels_1(i), trial_per_level*length(factor2_levels_1));
            conditions = [conditions; [factor1_levels_2(:) factor2_levels_3(:)]];
        end
    end
    save('.\temp\conditions.mat', 'conditions');
    exp_num = 1;
    start_trial = 1;
    save('.\temp\exp.mat', 'exp_num');
    save('.\temp\start_trial.mat', 'start_trial');
end

load('.\temp\conditions.mat')
load('.\temp\start_trial.mat');
load('.\temp\exp.mat');

if(exp_num == 1)
    for trial = start_trial :length(conditions)
        start_trial = trial;
        save('.\temp\start_trial.mat', 'start_trial');
        clc;
        
        factor1_level = conditions(trial, 1);
        factor2_level = conditions(trial, 2);
        if(factor2_level>0)
            F3_f_start = find_F3_f_start(F3_min_trans_f, F3_max_trans_f, factor2_level);
            F3_f_end  = F3_f_stable;
        else
            F3_f_start = -1;
            F3_f_end  = -1;
        end
        if(debug), disp([factor1_level factor2_level]); end
        
        TONE = zeros( tone_len, 2);
        
        switch factor1_level
            case 1
                [~, F3_trans_tone, ~] = gen_complex_tone(F3_f_start, F3_f_end, -1, trans_dur, -1, FS);
                TONE(transition_indx, right_ear) = F3_trans_tone;
            case 2
                [F1_cmplx_tone, ~, ~] = gen_complex_tone(F1_f_start, F1_f_end, F1_f_stable, trans_dur, stable_dur, FS);
                [F2_cmplx_tone, ~, ~] = gen_complex_tone(F2_f_start, F2_f_end, F2_f_stable, trans_dur, stable_dur, FS);
                [~, ~, F3_stable_tone] = gen_complex_tone(-1, -1, F3_f_stable, -1, stable_dur, FS);
                TONE(:, left_ear) = (F1_cmplx_tone + F2_cmplx_tone)/1;
                TONE(stable_indx, left_ear) = TONE(stable_indx, left_ear) + F3_stable_tone(:);
            case 3
                [~, F3_trans_tone, F3_stable_tone] = gen_complex_tone(F3_f_start, F3_f_end, F3_f_stable, trans_dur, stable_dur, FS);
                [F1_cmplx_tone, ~, ~] = gen_complex_tone(F1_f_start, F1_f_end, F1_f_stable, trans_dur, stable_dur, FS);
                [F2_cmplx_tone, ~, ~] = gen_complex_tone(F2_f_start, F2_f_end, F2_f_stable, trans_dur, stable_dur, FS);
                TONE(:, left_ear) = (F1_cmplx_tone + F2_cmplx_tone)/1;
                TONE(stable_indx, left_ear) = TONE(stable_indx, left_ear) + F3_stable_tone(:);
                TONE(transition_indx, right_ear) = TONE(transition_indx, right_ear) + F3_trans_tone(:);
            case 4
                [F1_cmplx_tone, ~, ~] = gen_complex_tone(F1_f_start, F1_f_end, F1_f_stable, trans_dur, stable_dur, FS);
                [F2_cmplx_tone, ~, ~] = gen_complex_tone(F2_f_start, F2_f_end, F2_f_stable, trans_dur, stable_dur, FS);
                [F3_cmplx_tone, ~, ~] = gen_complex_tone(F3_f_start, F3_f_end, F3_f_stable, trans_dur, stable_dur, FS);
                TONE = repmat((F1_cmplx_tone + F2_cmplx_tone + F3_cmplx_tone)/1, 1, 2);
        end
        
        stimuli_to_play = [silence; TONE; silence].*w;
        soundsc(stimuli_to_play, FS)
        pause(length(stimuli_to_play)/FS);
        if(debug)
            figure(1);
            spectrogram(stimuli_to_play(:,right_ear)+rand(size(stimuli_to_play(:,right_ear)))/10, 128, 32, 1024, 'yaxis');
            title('Right ear')
            set(gcf, 'WindowStyle', 'docked')
            figure(2);
            spectrogram(stimuli_to_play(:,left_ear)+rand(size(stimuli_to_play(:,left_ear)))/10, 128, 32, 1024, 'yaxis');
            title('Left ear')
            set(gcf, 'WindowStyle', 'docked')
        end
        
        valid_answer = 0;
        while(~valid_answer)
            disp('Which one did you hear?');
            disp('(1) /da/'); disp('(2) /ga/'); disp('(3) Chirp1'); disp('(4) Chirp2');
            disp('(5) /da/+Chirp1'); disp('(6) /da/+Chirp2'); disp('(7) /ga/+Chirp1'); disp('(8) /ga/+Chirp2')
            R = input('Enter your answer: ', 's');
            if(isempty(R) || length(R)>1)
                disp('invalid answer'); pause(1);
            else
                if( strcmpi(R, '1') || strcmpi(R, '2') || strcmpi(R, '3') || strcmpi(R, '4')  || strcmpi(R, '5') || strcmpi(R, '6') || strcmpi(R, '7') || strcmpi(R, '8')), valid_answer = 1; else; disp('invalid answer'); pause(1); end
            end
        end
        
        load('response_experiment1.mat');
        response = [response; {listener_code factor1_level factor2_level F3_f_start F3_f_end str2num(R)}];
        save('response_experiment1.mat', 'response');
        
    end
    
    exp_num  =2;
    save('.\temp\exp.mat', 'exp_num');
    start_trial = 1;
    save('.\temp\start_trial.mat', 'start_trial');
    recovery_mode = 0;
end
%% experiment 2
if(exp_num == 2)
    
    clc
    input('Press enter to begin experiment 2');
    
    factor1_levels = 1:6;
    factor2_levels = 1:2;
    factor1_levels_1 = factor1_levels(randperm(length(factor1_levels)));
    
    if(~recovery_mode)
        if(rand(1,1)<=.5)
            chirp_ear = left_ear;
            speech_ear = right_ear;
        else
            chirp_ear = right_ear;
            speech_ear = left_ear;
        end
        save('.\temp\chirp_ear.mat', 'chirp_ear');
        save('.\temp\speech_ear.mat', 'speech_ear');
    else
        load('.\temp\chirp_ear.mat');
        load('.\temp\speech_ear.mat');
    end
    
    if(~recovery_mode)
        conditions = [];
        for i = 1 : length(factor1_levels)
            
            factor2_levels_1 = factor2_levels(randperm(length(factor2_levels)));
            factor2_levels_2 = repelem(factor2_levels_1, trial_per_level);
            factor2_levels_3 = factor2_levels_2( randperm(length(factor2_levels_2)) );
            factor1_levels_2 = repelem(factor1_levels_1(i), trial_per_level*length(factor2_levels_1));
            conditions = [conditions; [factor1_levels_2(:) factor2_levels_3(:)]];
            
        end
        conditions = conditions(randperm(length(conditions)), :);
        save('.\temp\conditions.mat', 'conditions');
        exp_num = 2;
        start_trial = 1;
        save('.\temp\exp.mat', 'exp_num');
        save('.\temp\start_trial.mat', 'start_trial');
    end
    
    load('.\temp\conditions.mat')
    load('.\temp\start_trial.mat');
    load('.\temp\exp.mat');
    
    
    for trial = start_trial : length(conditions)
        start_trial = trial;
        save('.\temp\start_trial.mat', 'start_trial');
        clc;
        
        factor1_level = conditions(trial, 1); % formant transition difference. 1-4, 2-5, ...
        factor2_level = conditions(trial, 2);  %   stimuli type: chirp or speech
        if(debug), disp([factor1_level factor2_level]); end
        
        if(rand(1,1)<=0.5)
            F3_f_start_1 = find_F3_f_start(F3_min_trans_f, F3_max_trans_f, factor1_level);
            F3_f_start_2 = find_F3_f_start(F3_min_trans_f, F3_max_trans_f, factor1_level+3);
        else
            F3_f_start_2 = find_F3_f_start(F3_min_trans_f, F3_max_trans_f, factor1_level);
            F3_f_start_1 = find_F3_f_start(F3_min_trans_f, F3_max_trans_f, factor1_level+3);
        end
        F3_f_end  = F3_f_stable;
        if(debug), disp([round(F3_f_start_1) round(F3_f_start_2)]); end
        
        switch factor2_level
            case 1
                % chirp
                
                [~, F3_trans_tone_1, ~] = gen_complex_tone(F3_f_start_1, F3_f_end, -1, trans_dur, -1, FS);
                [~, F3_trans_tone_2, ~] = gen_complex_tone(F3_f_start_2, F3_f_end, -1, trans_dur, -1, FS);
                
                silence_between = zeros(0.05*FS, 1);
                TONE = zeros( 2*trans_dur*FS/1000 + length(silence_between), 2);
                TONE(:, chirp_ear) = [F3_trans_tone_1(:); silence_between; F3_trans_tone_2(:)];
            case 2
                % speech
                [F1_cmplx_tone, ~, ~] = gen_complex_tone(F1_f_start, F1_f_end, F1_f_stable, trans_dur, stable_dur, FS);
                [F2_cmplx_tone, ~, ~] = gen_complex_tone(F2_f_start, F2_f_end, F2_f_stable, trans_dur, stable_dur, FS);
                
                [F3_cmplx_tone_1, ~, ~] = gen_complex_tone(F3_f_start_1, F3_f_end, F3_f_stable, trans_dur, stable_dur, FS);
                [F3_cmplx_tone_2, ~, ~] = gen_complex_tone(F3_f_start_2, F3_f_end, F3_f_stable, trans_dur, stable_dur, FS);
                
                silence_between = zeros(0.25*FS, 1);
                TONE = zeros( 2*tone_len + length(silence_between), 2);
                TONE(:, speech_ear) = [(F1_cmplx_tone + F2_cmplx_tone + F3_cmplx_tone_1)/1; ...
                    silence_between; ...
                    (F1_cmplx_tone + F2_cmplx_tone + F3_cmplx_tone_2)/1];
                
        end
        
        silence = zeros(0.5*FS,2);
        w = hanning(size(TONE,1) + 2*length(silence));
        
        stimuli_to_play = [silence; TONE; silence].*w;
        soundsc(stimuli_to_play, FS)
        pause(length(stimuli_to_play)/FS);
        if(debug)
            figure(1);
            spectrogram(stimuli_to_play(:,right_ear)+rand(size(stimuli_to_play(:,right_ear)))/10, 128, 32, 1024, 'yaxis');
            title('Right ear')
            set(gcf, 'WindowStyle', 'docked')
            figure(2);
            spectrogram(stimuli_to_play(:,left_ear)+rand(size(stimuli_to_play(:,left_ear)))/10, 128, 32, 1024, 'yaxis');
            title('Left ear')
            set(gcf, 'WindowStyle', 'docked')
        end
        
        valid_answer = 0;
        while(~valid_answer)
            disp('Did the stimuli sound different');
            disp('(1) Yes'); disp('(2) No');
            R = input('Enter your answer: ', 's');
            if(isempty(R) || length(R)>1)
                disp('invalid answer'); pause(1);
            else
                if( strcmpi(R, '1') || strcmpi(R, '2')), valid_answer = 1; else; disp('invalid answer'); pause(1); end
            end
        end
        
        load('response_experiment2.mat');
        if(chirp_ear == left_ear)
            response = [response; {listener_code factor1_level factor2_level F3_f_start_1 F3_f_start_2 F3_f_end 'left' 'right' str2num(R)}];
        else
            response = [response; {listener_code factor1_level factor2_level F3_f_start_1 F3_f_start_2 F3_f_end 'right' 'left' str2num(R)}];
        end
        save('response_experiment2.mat', 'response');
        
    end
    
end
disp('Experiment is completed.')

clear