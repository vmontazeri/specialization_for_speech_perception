clear
close all
clc;

response = {'listener', 'factor1_level', 'factor2_level', 'F3_start', 'F3_end', 'response', 'correct_response'};
save('response_experiment1.mat', 'response');

response = {'listener', 'factor1_level', 'factor2_level', 'F3_f_start_1', 'F3_f_start_2', 'F3_f_end', 'chirp_ear', 'speech_ear', 'response', 'correct_response'};
save('response_experiment2.mat', 'response');

training_response = {'listener', 'factor1_level', 'factor2_level', 'F3_start', 'F3_end', 'response', 'correct_response'};
save('training_response.mat', 'training_response');

disp('Setup completed!');