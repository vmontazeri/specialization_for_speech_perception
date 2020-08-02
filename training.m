function training( debug, listener_code )

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

silence = zeros(0.5*FS,1);
tone_len = (trans_dur + stable_dur)*FS/1000;
w = hanning(tone_len + 2*length(silence));

transition_indx = 1:trans_dur*FS/1000;
stable_indx = trans_dur*FS/1000+1: tone_len;

trial_per_level = 5;
factor1_levels = 4;
factor2_levels = [1:4];
% factor1_levels_1 = factor1_levels(randperm(length(factor1_levels)));
factor1_levels_1 = factor1_levels;

conditions = [];
for i = 1 : length(factor1_levels)
    if(factor1_levels(i)==2)        
        factor1_levels_2 = repelem(factor1_levels_1(i), trial_per_level);
        conditions = [conditions; [factor1_levels_2(:) -1*ones(trial_per_level,1)]];
    else
        factor2_levels_1 = factor2_levels(randperm(length(factor2_levels)));
        factor2_levels_2 = repelem(factor2_levels_1, trial_per_level);
%         factor2_levels_3 = factor2_levels_2( randperm(length(factor2_levels_2)) );
        factor2_levels_3 = factor2_levels_2;
        factor1_levels_2 = repelem(factor1_levels_1(i), trial_per_level*length(factor2_levels_1));
        conditions = [conditions; [factor1_levels_2(:) factor2_levels_3(:)]];
    end
end
conditions = [conditions; conditions];

clc;
disp('Before doing the experiment, you are going to have a brief training session to become familiar with the stimuli.')
disp('Please note that during the traning, you''ll recevie feedback once you submit your answer!');
input('Press enter to proceed');

for trial = 1 : length(conditions)
    clc;
    
    factor1_level = conditions(trial, 1);
    factor2_level = conditions(trial, 2);
    
    if(debug), disp([factor1_level factor2_level]); end
    
    TONE = zeros( tone_len, 1);
    
    switch factor2_level
        case 1 % da
             F3_f_start = find_F3_f_start(F3_min_trans_f, F3_max_trans_f, 9);
             F3_f_end  = F3_f_stable;
            [F1_cmplx_tone, ~, ~] = gen_complex_tone(F1_f_start, F1_f_end, F1_f_stable, trans_dur, stable_dur, FS);
            [F2_cmplx_tone, ~, ~] = gen_complex_tone(F2_f_start, F2_f_end, F2_f_stable, trans_dur, stable_dur, FS);
            [F3_cmplx_tone, ~, ~] = gen_complex_tone(F3_f_start, F3_f_end, F3_f_stable, trans_dur, stable_dur, FS);            
            TONE = (F1_cmplx_tone + F2_cmplx_tone + F3_cmplx_tone)/1;               
        case 2 % ga
            F3_f_start = find_F3_f_start(F3_min_trans_f, F3_max_trans_f, 1);
            F3_f_end  = F3_f_stable;
            [F1_cmplx_tone, ~, ~] = gen_complex_tone(F1_f_start, F1_f_end, F1_f_stable, trans_dur, stable_dur, FS);
            [F2_cmplx_tone, ~, ~] = gen_complex_tone(F2_f_start, F2_f_end, F2_f_stable, trans_dur, stable_dur, FS);
            [F3_cmplx_tone, ~, ~] = gen_complex_tone(F3_f_start, F3_f_end, F3_f_stable, trans_dur, stable_dur, FS);            
            TONE = (F1_cmplx_tone + F2_cmplx_tone + F3_cmplx_tone)/1;   
        case 3 %chirp1
            F3_f_start = find_F3_f_start(F3_min_trans_f, F3_max_trans_f, 1);
            F3_f_end  = F3_f_stable;            
            [~, F3_trans_tone, ~] = gen_complex_tone(F3_f_start, F3_f_end, F3_f_stable, trans_dur, stable_dur, FS);            
            TONE(transition_indx, 1) = F3_trans_tone(:); 
        case 4 %chirp2
            F3_f_start = find_F3_f_start(F3_min_trans_f, F3_max_trans_f, 9);
            F3_f_end  = F3_f_stable;            
            [~, F3_trans_tone, ~] = gen_complex_tone(F3_f_start, F3_f_end, F3_f_stable, trans_dur, stable_dur, FS);            
            TONE(transition_indx, 1) = F3_trans_tone(:); 
    end
    
    correct_answer = factor2_level;
    
    stimuli_to_play = [silence; TONE; silence].*w;
    soundsc(stimuli_to_play, FS)
    pause(length(stimuli_to_play)/FS);
    if(debug)
        figure(1); 
        spectrogram(stimuli_to_play(:,1)+rand(size(stimuli_to_play(:,1)))/10, 128, 32, 1024, 'yaxis'); 
        title('Right ear')
        set(gcf, 'WindowStyle', 'docked')
        figure(2); 
        spectrogram(stimuli_to_play(:,1)+rand(size(stimuli_to_play(:,1)))/10, 128, 32, 1024, 'yaxis'); 
        title('Left ear')
        set(gcf, 'WindowStyle', 'docked')
    end
    
    valid_answer = 0;
    while(~valid_answer)
        disp('Which one did you hear?'); 
        disp('(1) /da/'); disp('(2) /ga/'); disp('(3) Chirp1'); disp('(4) Chirp2');         
        R = input('Enter your answer: ', 's');
        if(length(R)>1 || isempty(R))
            disp('invalid answer'); pause(1);
        else
            if( strcmpi(R, '1') || strcmpi(R, '2') || strcmpi(R, '3') || strcmpi(R, '4')), valid_answer = 1; 
            else; disp('invalid answer'); pause(1); 
            end
        end
    end
    
    if(str2num(R) == correct_answer)
        disp('CORRECT :)');
    else
        switch correct_answer
            case 1
                correct_answer_details = '/da/';
            case 2
                correct_answer_details = '/ga/';
            case 3
                correct_answer_details = 'Chirp 1';
            case 4
                correct_answer_details = 'Chirp 2';
        end                
        disp( ['INCORRECT :( Correct answer was ' correct_answer_details] );
    end
    
    disp('')
    input('Press enter to proceed')
    
    load('training_response.mat');
    training_response = [training_response; {listener_code factor1_level factor2_level F3_f_start F3_f_end str2num(R)}];
    save('training_response.mat', 'training_response');

end
