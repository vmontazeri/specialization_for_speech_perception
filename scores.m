clear
close all
clc

load('response_experiment1.mat');

res = zeros(length(response)-1, 8);
total = zeros(8, 1);
bias = total;

for i = 2 : length(response)
    
    row = [0 cell2mat(response(i,2:end))];
    part_name = char( response(i,1) );
            
    bias(row(end)) = bias(row(end)) + 1;
      
    if(row(2) == 1 && row(3)<5), corr = 3; total(corr) = total(corr)+1; end
    if(row(2) == 1 && row(3)>5), corr = 4; total(corr) = total(corr)+1; end
    
    if(row(2) == 2), corr = 1;  total(corr) = total(corr)+1; end    
    
    if(row(2) == 3 && row(3)<5), corr = 7; total(corr) = total(corr)+1; end
    if(row(2) == 3 && row(3)>5), corr = 6; total(corr) = total(corr)+1; end
    
    if(row(2) == 4 && row(3)<5), corr = 2; total(corr) = total(corr)+1; end
    if(row(2) == 4 && row(3)>5), corr = 1; total(corr) = total(corr)+1; end
    
    if(row(end) == corr), res(i,corr) = 1; end
        
end

disp('Bias:')
100*(bias / sum(bias))'


100*sum(res) ./ total'
%%
clear
close all
clc

load('response_experiment2.mat');

res = zeros(2, 6);
total = zeros(4, 1);
bias = total;

for i = 2 : length(response)
    
    row = cell2mat(response(i,[2 3 9]));
    part_name = char( response(i,1) );
            
    if(contains(part_name, 'MKFNH1'))        
        disp(part_name);
        
        one_flag = 0; two_flag = 0;
        if(row(end) == 1) 
            one_flag=1; 
        end
        if(row(end) == 2) 
            two_flag=1; 
        end
        
        if(one_flag), row(end) = 2; end
        if(two_flag), row(end) = 1; end        
        
    end


%     bias(row(end)) = bias(row(end)) + 1;
      
    if(row(2) == 1), corr = 1; total(corr) = total(corr)+1; end   
    
    if(row(2) == 2), corr = 2;  total(corr) = total(corr)+1; end    
    
    if(row(end) == 1), res(corr,row(1)) = res(corr,row(1))+1; end
        
end

plot(1:6, 100*res'/60, '-s', 'LineWidth', 2)
legend([{'non-speech'}, {'speech'}])
xticks(0:7)
xticklabels({'', '1-4','2-5','3-6','4-7','5-8','6-9', ''})
xlim([0 7])
grid on
ylabel('Percent correct')
xlabel('Stimulus pair')