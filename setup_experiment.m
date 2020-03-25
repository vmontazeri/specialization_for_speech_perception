clear
close all
clc;

response = {'listener', 'factor1_level', 'factor2_level', 'F3_start', 'F3_end', 'response'};
save('response.mat', 'response');

training_response = {'listener', 'factor1_level', 'factor2_level', 'F3_start', 'F3_end', 'response'};
save('training_response.mat', 'training_response');

disp('Setup completed!');