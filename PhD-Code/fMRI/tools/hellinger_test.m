
interval = -1500:1500;
mn = 0;
stdeviation1A = 200;
stdeviation2A = 250;
stdeviation1B = 400;
stdeviation2B = 450;

Y1A = normpdf(interval,mn,stdeviation1A);
Y2A = normpdf(interval,mn,stdeviation2A);
Y1B = normpdf(interval,mn,stdeviation1B);
Y2B = normpdf(interval,mn,stdeviation2B);

max_y = max([Y1A, Y2A, Y1B, Y2B]);

x = 500;
y = 0.9 * max_y;

figure;

subplot(1,3,1);
plot(interval,Y1A,'b');
hold on
plot(interval,Y2A,'r');
ylim([0 max_y]);
text(x,y,strcat('Hellinger:',num2str(sqrt(1 - sum( sqrt(Y1A.*Y2A) )))));

subplot(1,3,2);
plot(interval,Y1B,'g');
hold on
plot(interval,Y2B,'y');
ylim([0 max_y]);
text(x,y,strcat('Hellinger:',num2str(sqrt(1 - sum( sqrt(Y1B.*Y2B) )))));

subplot(1,3,3);
plot(interval,Y1A,'b');
hold on
plot(interval,Y1B,'g');
ylim([0 max_y]);
text(x,y,strcat('Hellinger:',num2str(sqrt(1 - sum( sqrt(Y1A.*Y1B) )))));
