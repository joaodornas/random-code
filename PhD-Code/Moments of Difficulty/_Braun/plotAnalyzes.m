
f = figure;

MarkerSize = 10;
LineWidth = 2;

nWindows = length(Ecolor);

subplot(3,1,1);

plot(Ecolor,'b');
hold on

plot(kECmax,Ecolor(kECmax),'ro','LineWidth',LineWidth,'MarkerSize',MarkerSize);
plot(kECmin,Ecolor(kECmin),'ro','LineWidth',LineWidth,'MarkerSize',MarkerSize);

xlim([0 nWindows]);
title('Color Entropy');
xlabel('windows');
ylabel('Entropy');

subplot(3,1,2);

%g = figure;

plot(Emotion,'b');
hold on

plot(kEMmax,Emotion(kEMmax),'ro','LineWidth',LineWidth,'MarkerSize',MarkerSize);
plot(kEMmin,Emotion(kEMmin),'ro','LineWidth',LineWidth,'MarkerSize',MarkerSize);

xlim([0 nWindows]);
title('Motion Entropy');
xlabel('windows');
ylabel('Entropy');

subplot(3,1,3);

%f = figure;

area = squeeze(areaComb(minArea.color,minArea.comb,:));

plot(area,'b');
hold on

idx_max = find(area == max(area));

plot(minArea.kf,area(minArea.kf),'ro','LineWidth',LineWidth,'MarkerSize',MarkerSize);
plot(idx_max,area(idx_max),'ro','LineWidth',LineWidth,'MarkerSize',MarkerSize);

xlim([0 nWindows]);
title('Area');
xlabel('windows');
ylabel('Area');

print(f,'-djpeg',strcat(datafilename(1:end-4),'-','Color-Motion-Area','.jpeg'));
print(f,'-dpdf',strcat(datafilename(1:end-4),'-','Color-Motion-Area','.pdf'));
print(f,'-depsc',strcat(datafilename(1:end-4),'-','Color-Motion-Area','.eps'));