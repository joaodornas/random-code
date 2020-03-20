clear all
close all

dio = digitalio('parallel','LPT1');
addline(dio,[1],'in');
g = [];
start = GetSecs;
i = 1;
n =[start];
figure
while ~KbCheck
%     while getvalue(dio)
        g = [g;getvalue(dio)];
        if length(g)>101
            plot(g(end-100:end,1))
            ylim([0;1.1]);
            drawnow
        end
%     end
%     n = [n,GetSecs];
%     while bin2dec(num2str(getvalue(dio))) > 1
%         g = [g;getvalue(dio),GetSecs]];
%     end
end

plot(g(:,1));
    
