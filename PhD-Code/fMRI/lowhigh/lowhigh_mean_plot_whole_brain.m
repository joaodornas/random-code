
function lowhigh_mean_plot_whole_brain

settings_jan_0502;

get_at_this_preprocessed_step = settings.folders.realigned.name;
prefix_for_the_preprocessed_step = settings.folders.realigned.prefix;

lowhigh_load_all_data;

k = 3; %% polynomial order

frame_size = 120 - 1;
doWithGolayAndPlot(MOT4Run1,MOT4Run2,MOT2Run1,MOT2Run2,RestingStateRun1,RestingStateRun2,k,frame_size,settings);

frame_size = 120/2 - 1;
doWithGolayAndPlot(MOT4Run1,MOT4Run2,MOT2Run1,MOT2Run2,RestingStateRun1,RestingStateRun2,k,frame_size,settings);

frame_size = 120/4 - 1;
doWithGolayAndPlot(MOT4Run1,MOT4Run2,MOT2Run1,MOT2Run2,RestingStateRun1,RestingStateRun2,k,frame_size,settings);

frame_size = 6 - 1;
doWithGolayAndPlot(MOT4Run1,MOT4Run2,MOT2Run1,MOT2Run2,RestingStateRun1,RestingStateRun2,k,frame_size,settings);

frame_size = 120 - 1;
doRemovingSlowDrifts(MOT4Run1,MOT4Run2,MOT2Run1,MOT2Run2,RestingStateRun1,RestingStateRun2,k,frame_size,settings);

doWithOutGolayAndPlot(MOT4Run1,MOT4Run2,MOT2Run1,MOT2Run2,RestingStateRun1,RestingStateRun2,settings);

doBOLDSignalAndPlot(MOT4Run1,MOT4Run2,MOT2Run1,MOT2Run2,RestingStateRun1,RestingStateRun2,settings);

return

function doWithGolayAndPlot(MOT4Run1,MOT4Run2,MOT2Run1,MOT2Run2,RestingStateRun1,RestingStateRun2,k,frame_size,settings)

disp('doWithGolayAndPlot');

fontSize = 6;

MOT4Run1sgolay = sgolayfilt(MOT4Run1,k,frame_size,[],4);
MOT4Run2sgolay = sgolayfilt(MOT4Run2,k,frame_size,[],4);

MOT4Run1ZStgolay = zscore(MOT4Run1sgolay,0,4);
MOT4Run2ZStgolay = zscore(MOT4Run2sgolay,0,4);

MOT4Run1ZStmeangolay = mean_4D(MOT4Run1ZStgolay);
MOT4Run2ZStmeangolay = mean_4D(MOT4Run2ZStgolay);

MOT2Run1sgolay = sgolayfilt(MOT2Run1,k,frame_size,[],4);
MOT2Run2sgolay = sgolayfilt(MOT2Run2,k,frame_size,[],4);

MOT2Run1ZStgolay = zscore(MOT2Run1sgolay,0,4);
MOT2Run2ZStgolay = zscore(MOT2Run2sgolay,0,4);

MOT2Run1ZStmeangolay = mean_4D(MOT2Run1ZStgolay);
MOT2Run2ZStmeangolay = mean_4D(MOT2Run2ZStgolay);

RestingStateRun1sgolay = sgolayfilt(RestingStateRun1,k,frame_size,[],4);
RestingStateRun2sgolay = sgolayfilt(RestingStateRun2,k,frame_size,[],4);

RestingStateRun1ZStgolay = zscore(RestingStateRun1sgolay,0,4);
RestingStateRun2ZStgolay = zscore(RestingStateRun2sgolay,0,4);

RestingStateRun1ZStmeangolay = mean_4D(RestingStateRun1ZStgolay);
RestingStateRun2ZStmeangolay = mean_4D(RestingStateRun2ZStgolay);

disp('Plot - mit golay');

h = figure;

subplot(3,2,1);
plot(MOT4Run1ZStmeangolay,'r');
hold on
plot(MOT4Run2ZStmeangolay,'b');
legend('Run1', 'Run2');
title('MOT4');
xlabel('TR');
ylabel(strcat('BOLD - Savitzky-Golay (',int2str(frame_size),') + z-scored'),'FontSize',fontSize);

hold on

subplot(3,2,3);
plot(MOT2Run1ZStmeangolay,'r');
hold on
plot(MOT2Run2ZStmeangolay,'b');
legend('Run1', 'Run2');
title('MOT2');
xlabel('TR');
ylabel(strcat('BOLD - Savitzky-Golay (',int2str(frame_size),') + z-scored'),'FontSize',fontSize);

hold on

subplot(3,2,5);
plot(RestingStateRun1ZStmeangolay,'r');
hold on
plot(RestingStateRun2ZStmeangolay,'b');
legend('Run1', 'Run2');
title('Resting State');
xlabel('TR');
ylabel(strcat('BOLD - Savitzky-Golay (',int2str(frame_size),') + z-scored'),'FontSize',fontSize);

hold on

subplot(3,2,2);
plot(MOT4Run1ZStmeangolay,'r');
hold on
plot(RestingStateRun1ZStmeangolay,'b');
legend('MOT4 - Run1', 'Resting State - Run1');
title('MOT4 x Resting State');
xlabel('TR');
ylabel(strcat('BOLD - Savitzky-Golay (',int2str(frame_size),') + z-scored'),'FontSize',fontSize);

hold on

subplot(3,2,4);
plot(MOT2Run1ZStmeangolay,'r');
hold on
plot(RestingStateRun1ZStmeangolay,'b');
legend('MOT2 - Run1', 'Resting State - Run1');
title('MOT2 x Resting State');
xlabel('TR');
ylabel(strcat('BOLD - Savitzky-Golay (',int2str(frame_size),') + z-scored'),'FontSize',fontSize);

hold on

subplot(3,2,6);
plot(MOT4Run1ZStmeangolay,'r');
hold on
plot(MOT2Run1ZStmeangolay,'b');
legend('MOT4 - Run1', 'MOT2 - Run1');
title('MOT4 x MOT2');
xlabel('TR');
ylabel(strcat('BOLD - Savitzky-Golay (',int2str(frame_size),') + z-scored'),'FontSize',fontSize);

print(h, '-djpeg', strcat(settings.folders.experiment,'-',settings.subject,'-','mean-whole-brain-z-scored-golay-',int2str(frame_size)));

return

function doRemovingSlowDrifts(MOT4Run1,MOT4Run2,MOT2Run1,MOT2Run2,RestingStateRun1,RestingStateRun2,k,frame_size,settings)

disp('doRemovingSlowDrifts');

fontSize = 6;

MOT4Run1sgolay = sgolayfilt(MOT4Run1,k,frame_size,[],4);
MOT4Run2sgolay = sgolayfilt(MOT4Run2,k,frame_size,[],4);

MOT4Run1ZStgolay = zscore(MOT4Run1sgolay,0,4);
MOT4Run2ZStgolay = zscore(MOT4Run2sgolay,0,4);

MOT4Run1ZStmeangolay = mean_4D(MOT4Run1ZStgolay);
MOT4Run2ZStmeangolay = mean_4D(MOT4Run2ZStgolay);

MOT4Run1ZSt = zscore(MOT4Run1,0,4);
MOT4Run2ZSt = zscore(MOT4Run2,0,4);

MOT4Run1ZStmean = mean_4D(MOT4Run1ZSt);
MOT4Run2ZStmean = mean_4D(MOT4Run2ZSt);

MOT4Run1removed = MOT4Run1ZStmean - MOT4Run1ZStmeangolay;
MOT4Run2removed = MOT4Run2ZStmean - MOT4Run2ZStmeangolay;

MOT2Run1sgolay = sgolayfilt(MOT2Run1,k,frame_size,[],4);
MOT2Run2sgolay = sgolayfilt(MOT2Run2,k,frame_size,[],4);

MOT2Run1ZStgolay = zscore(MOT2Run1sgolay,0,4);
MOT2Run2ZStgolay = zscore(MOT2Run2sgolay,0,4);

MOT2Run1ZStmeangolay = mean_4D(MOT2Run1ZStgolay);
MOT2Run2ZStmeangolay = mean_4D(MOT2Run2ZStgolay);

MOT2Run1ZSt = zscore(MOT2Run1,0,4);
MOT2Run2ZSt = zscore(MOT2Run2,0,4);

MOT2Run1ZStmean = mean_4D(MOT2Run1ZSt);
MOT2Run2ZStmean = mean_4D(MOT2Run2ZSt);

MOT2Run1removed = MOT2Run1ZStmean - MOT2Run1ZStmeangolay;
MOT2Run2removed = MOT2Run2ZStmean - MOT2Run2ZStmeangolay;

RestingStateRun1sgolay = sgolayfilt(RestingStateRun1,k,frame_size,[],4);
RestingStateRun2sgolay = sgolayfilt(RestingStateRun2,k,frame_size,[],4);

RestingStateRun1ZStgolay = zscore(RestingStateRun1sgolay,0,4);
RestingStateRun2ZStgolay = zscore(RestingStateRun2sgolay,0,4);

RestingStateRun1ZStmeangolay = mean_4D(RestingStateRun1ZStgolay);
RestingStateRun2ZStmeangolay = mean_4D(RestingStateRun2ZStgolay);

RestingStateRun1ZSt = zscore(RestingStateRun1,0,4);
RestingStateRun2ZSt = zscore(RestingStateRun2,0,4);

RestingStateRun1ZStmean = mean_4D(RestingStateRun1ZSt);
RestingStateRun2ZStmean = mean_4D(RestingStateRun2ZSt);

RestingStateRun1removed = RestingStateRun1ZStmean - RestingStateRun1ZStmeangolay;
RestingStateRun2removed = RestingStateRun2ZStmean - RestingStateRun2ZStmeangolay;

disp('Plot - removed low-freq drifts');

h = figure;

subplot(3,2,1);
plot(MOT4Run1removed,'r');
hold on
plot(MOT4Run2removed,'b');
legend('Run1', 'Run2');
title('MOT4');
xlabel('TR');
ylabel(strcat('BOLD - removed low-freq drifts + z-scored'),'FontSize',fontSize);

hold on

subplot(3,2,3);
plot(MOT2Run1removed,'r');
hold on
plot(MOT2Run2removed,'b');
legend('Run1', 'Run2');
title('MOT2');
xlabel('TR');
ylabel(strcat('BOLD - removed low-freq drifts + z-scored'),'FontSize',fontSize);

hold on

subplot(3,2,5);
plot(RestingStateRun1removed,'r');
hold on
plot(RestingStateRun2removed,'b');
legend('Run1', 'Run2');
title('Resting State');
xlabel('TR');
ylabel(strcat('BOLD - removed low-freq drifts + z-scored'),'FontSize',fontSize);

hold on

subplot(3,2,2);
plot(MOT4Run1removed,'r');
hold on
plot(RestingStateRun1removed,'b');
legend('MOT4 - Run1', 'Resting State - Run1');
title('MOT4 x Resting State');
xlabel('TR');
ylabel(strcat('BOLD - removed low-freq drifts + z-scored'),'FontSize',fontSize);

hold on

subplot(3,2,4);
plot(MOT2Run1removed,'r');
hold on
plot(RestingStateRun1removed,'b');
legend('MOT2 - Run1', 'Resting State - Run1');
title('MOT2 x Resting State');
xlabel('TR');
ylabel(strcat('BOLD - removed low-freq drifts + z-scored'),'FontSize',fontSize);

hold on

subplot(3,2,6);
plot(MOT4Run1removed,'r');
hold on
plot(MOT2Run1removed,'b');
legend('MOT4 - Run1', 'MOT2 - Run1');
title('MOT4 x MOT2');
xlabel('TR');
ylabel(strcat('BOLD - removed low-freq drifts + z-scored'),'FontSize',fontSize);

print(h, '-djpeg', strcat(settings.folders.experiment,'-',settings.subject,'-','mean-whole-brain-z-scored-removed-low-freq-drifts'));

return

function doWithOutGolayAndPlot(MOT4Run1,MOT4Run2,MOT2Run1,MOT2Run2,RestingStateRun1,RestingStateRun2,settings)

disp('doWithOutGolayAndPlot');

fontSize = 6;

MOT4Run1ZSt = zscore(MOT4Run1,0,4);
MOT4Run2ZSt = zscore(MOT4Run2,0,4);

MOT4Run1ZStmean = mean_4D(MOT4Run1ZSt);
MOT4Run2ZStmean = mean_4D(MOT4Run2ZSt);

MOT2Run1ZSt = zscore(MOT2Run1,0,4);
MOT2Run2ZSt = zscore(MOT2Run2,0,4);

MOT2Run1ZStmean = mean_4D(MOT2Run1ZSt);
MOT2Run2ZStmean = mean_4D(MOT2Run2ZSt);

RestingStateRun1ZSt = zscore(RestingStateRun1,0,4);
RestingStateRun2ZSt = zscore(RestingStateRun2,0,4);

RestingStateRun1ZStmean = mean_4D(RestingStateRun1ZSt);
RestingStateRun2ZStmean = mean_4D(RestingStateRun2ZSt);

disp('Plot - ohne golay');

h = figure;

subplot(3,2,1);
plot(MOT4Run1ZStmean,'r');
hold on
plot(MOT4Run2ZStmean,'b');
legend('Run1', 'Run2');
title('MOT4');
xlabel('TR');
ylabel('BOLD - z-scored','FontSize',fontSize);

hold on

subplot(3,2,3);
plot(MOT2Run1ZStmean,'r');
hold on
plot(MOT2Run2ZStmean,'b');
legend('Run1', 'Run2');
title('MOT2');
xlabel('TR');
ylabel('BOLD - z-scored','FontSize',fontSize);

hold on

subplot(3,2,5);
plot(RestingStateRun1ZStmean,'r');
hold on
plot(RestingStateRun2ZStmean,'b');
legend('Run1', 'Run2');
title('Resting State');
xlabel('TR');
ylabel('BOLD - z-scored','FontSize',fontSize);

hold on

subplot(3,2,2);
plot(MOT4Run1ZStmean,'r');
hold on
plot(RestingStateRun1ZStmean,'b');
legend('MOT4 - Run1', 'Resting State - Run1');
title('MOT4 x Resting State');
xlabel('TR');
ylabel('BOLD - z-scored','FontSize',fontSize);

hold on

subplot(3,2,4);
plot(MOT2Run1ZStmean,'r');
hold on
plot(RestingStateRun1ZStmean,'b');
legend('MOT2 - Run1', 'Resting State - Run1');
title('MOT2 x Resting State');
xlabel('TR');
ylabel('BOLD - z-scored','FontSize',fontSize);

hold on

subplot(3,2,6);
plot(MOT4Run1ZStmean,'r');
hold on
plot(MOT2Run1ZStmean,'b');
legend('MOT4 - Run1', 'MOT2 - Run1');
title('MOT4 x MOT2');
xlabel('TR');
ylabel('BOLD - z-scored','FontSize',fontSize);

print(h, '-djpeg', strcat(settings.folders.experiment,'-',settings.subject,'-','mean-whole-brain-z-scored'));

return

function doBOLDSignalAndPlot(MOT4Run1,MOT4Run2,MOT2Run1,MOT2Run2,RestingStateRun1,RestingStateRun2,settings)

disp('doBOLDSignalAndPlot');

fontSize = 6;

MOT4Run1mean = mean_4D(MOT4Run1);
MOT4Run2mean = mean_4D(MOT4Run2);

MOT2Run1mean = mean_4D(MOT2Run1);
MOT2Run2mean = mean_4D(MOT2Run2);

RestingStateRun1mean = mean_4D(RestingStateRun1);
RestingStateRun2mean = mean_4D(RestingStateRun2);

disp('Plot - BOLD Signal');

h = figure;

subplot(3,2,1);
plot(MOT4Run1mean,'r');
hold on
plot(MOT4Run2mean,'b');
legend('Run1', 'Run2');
title('MOT4');
xlabel('TR');
ylabel('BOLD signal','FontSize',fontSize);

hold on

subplot(3,2,3);
plot(MOT2Run1mean,'r');
hold on
plot(MOT2Run2mean,'b');
legend('Run1', 'Run2');
title('MOT2');
xlabel('TR');
ylabel('BOLD signal','FontSize',fontSize);

hold on

subplot(3,2,5);
plot(RestingStateRun1mean,'r');
hold on
plot(RestingStateRun2mean,'b');
legend('Run1', 'Run2');
title('Resting State');
xlabel('TR');
ylabel('BOLD signal','FontSize',fontSize);

hold on

subplot(3,2,2);
plot(MOT4Run1mean,'r');
hold on
plot(RestingStateRun1mean,'b');
legend('MOT4 - Run1', 'Resting State - Run1');
title('MOT4 x Resting State');
xlabel('TR');
ylabel('BOLD signal','FontSize',fontSize);

hold on

subplot(3,2,4);
plot(MOT2Run1mean,'r');
hold on
plot(RestingStateRun1mean,'b');
legend('MOT2 - Run1', 'Resting State - Run1');
title('MOT2 x Resting State');
xlabel('TR');
ylabel('BOLD signal','FontSize',fontSize);

hold on

subplot(3,2,6);
plot(MOT4Run1mean,'r');
hold on
plot(MOT2Run1mean,'b');
legend('MOT4 - Run1', 'MOT2 - Run1');
title('MOT4 x MOT2');
xlabel('TR');
ylabel('BOLD signal','FontSize',fontSize);

print(h, '-djpeg', strcat(settings.folders.experiment,'-',settings.subject,'-','mean-whole-brain-BOLD-signal'));

return



