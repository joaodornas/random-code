
for t=1:1000
H1(t)=normrnd(0,0.05);
H2(t)=normrnd(0,0.10); 
H3(t)=normrnd(0,0.30);
end

map = brewermap(3,'Set1'); 
figure
histf(H1,-1.3:.01:1.3,'facecolor',map(1,:),'facealpha',.5,'edgecolor','none')
hold on
histf(H2,-1.3:.01:1.3,'facecolor',map(2,:),'facealpha',.5,'edgecolor','none')
histf(H3,-1.3:.01:1.3,'facecolor',map(3,:),'facealpha',.5,'edgecolor','none')
box off
axis tight
legalpha('H1','H2','H3','location','northwest')
legend boxoff

