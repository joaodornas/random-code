%function normalizePSTHforRasterPlot(nConditions,nTrialsPerCondition)

% normalizes PSTH generated with m_convolute for plotting with raster
% inputs are number of conditions and number of trials per condition.
% normalizePSTHforRasterPlot(nConditions,nTrialsPerCondition)

max_resp=zeros(nConditions,1);
for i=1:nConditions
    psth=eval(strcat('PSTH',int2str(i)));
    max_resp(i)=max(psth);
end
the_max=max(max_resp);
for i=1:nConditions
    psth=eval(strcat('PSTH',int2str(i)));
    psth_norm=(psth./the_max).*nTrialsPerCondition;
    assignin('base',strcat('psth_norm',int2str(i)),psth_norm);
end