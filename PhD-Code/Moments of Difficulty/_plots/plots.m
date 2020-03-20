eventsBalls;

figure

subplot(3,1,1);
plot(1:nFrames,vi(selectedBalls(1),1:nFrames));
title('Speed Ball 1');
ylim([0 maxvi]);

subplot(3,1,2);
plot(1:nFrames,vi(selectedBalls(2),1:nFrames));
title('Speed Ball 2');
ylim([0 maxvi]);

subplot(3,1,3);
plot(1:nFrames,vi(selectedBalls(3),1:nFrames));
title('Speed Ball 3');
ylim([0 maxvi]);

hold on

figure

subplot(3,1,1);
plot(1:nFrames,oi(selectedBalls(1),1:nFrames));
title('Angular Speed Ball 1');
ylim([0 maxoi]);

subplot(3,1,2);
plot(1:nFrames,oi(selectedBalls(2),1:nFrames));
title('Angular Speed Ball 2');
ylim([0 maxoi]);

subplot(3,1,3);
plot(1:nFrames,oi(selectedBalls(3),1:nFrames));
title('Angular Speed Ball 3');
ylim([0 maxoi]);

hold on

figure

plot(1:size(areaBalls,2),areaBalls(areaComb,:));
title('Area');

hold on

figure

plot(1:length(eventsArea), eventsArea,'bo');
title('Events Area');

hold on

figure

subplot(3,1,1);
plot(1:size(eventsBallDisplacement,2), eventsBallDisplacement(1,:),'bo');
title('Events Displacement Ball 1');

subplot(3,1,2);
plot(1:size(eventsBallDisplacement,2), eventsBallDisplacement(2,:),'bo');
title('Events Displacement Ball 2');

subplot(3,1,3);
plot(1:size(eventsBallDisplacement,2), eventsBallDisplacement(3,:),'bo');
title('Events Displacement Ball 3');

hold on

figure

subplot(3,1,1);
plot(1:size(eventsBallColorMatch,2), eventsBallColorMatch(1,:),'bo');
title('Events Color Match Ball 1');

subplot(3,1,2);
plot(1:size(eventsBallColorMatch,2), eventsBallColorMatch(2,:),'bo');
title('Events Color Match Ball 2');

subplot(3,1,3);
plot(1:size(eventsBallColorMatch,2), eventsBallColorMatch(3,:),'bo');
title('Events Color Match Ball 3');

hold on

figure

subplot(3,1,1);
plot(1:size(eventsBallMixing,3), sum(squeeze(eventsBallMixing(selectedBalls(1),:,:))),'bo');
title('Events Mixing Ball 1');

subplot(3,1,2);
plot(1:size(eventsBallMixing,3), sum(squeeze(eventsBallMixing(selectedBalls(2),:,:))),'bo');
title('Events Mixing Ball 2');

subplot(3,1,3);
plot(1:size(eventsBallMixing,3), sum(squeeze(eventsBallMixing(selectedBalls(3),:,:))),'bo');
title('Events Mixing Ball 3');




