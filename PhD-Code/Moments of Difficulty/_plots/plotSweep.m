ball = input('ball: ', 's');

DataFile= input('Please enter the Parameters file name: ', 's');

load(DataFile);

figure

for i=1:16
    
    subplot(8,2,i)

    plot(1:size(sweep,3),squeeze(sweep(str2num(ball),i,:)));
    title(strcat('Sweep per Ball:',int2str(i)));
    ylim([0 max(sweep(str2num(ball),i,:))]);
    
end




