function plotPM = plotPM(stdM,meanM,line)


PM = meanM + stdM;
MM = meanM - stdM;

len = size(PM,2);
x = 1:len;
y1 = PM(line,:);
y2 = MM(line,:);

plot(x,y1,x,y2)


end