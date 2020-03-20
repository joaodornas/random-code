figure

subplot(3,1,1);
plot(1:nFrames,vi(1,1:nFrames));
title('Speed Ball 1');
ylim([0 maxvi]);

subplot(3,1,2);
plot(1:nFrames,vi(2,1:nFrames));
title('Speed Ball 2');
ylim([0 maxvi]);

subplot(3,1,3);
plot(1:nFrames,vi(3,1:nFrames));
title('Speed Ball 3');
ylim([0 maxvi]);

hold on

figure

subplot(3,1,1);
plot(1:nFrames,oi(selectedBalls(1),1:nFrames));
title('Angular Speed Ball 1');

subplot(3,1,2);
plot(1:nFrames,oi(selectedBalls(2),1:nFrames));
title('Angular Speed Ball 2');

subplot(3,1,3);
plot(1:nFrames,oi(selectedBalls(3),1:nFrames));
title('Angular Speed Ball 3');

hold on

figure

plot(1:nFrames,areaBalls);
title('Area');

hold on

figure

subplot(3,1,1);
plot(1:nFrames,ballColorMatch(1,1:nFrames),'bo');
title('Color Match Ball 1');

subplot(3,1,2);
plot(1:nFrames,ballColorMatch(2,1:nFrames),'bo');
title('Color Match Ball 2');

subplot(3,1,3);
plot(1:nFrames,ballColorMatch(3,1:nFrames),'bo');
title('Color Match Ball 3');

