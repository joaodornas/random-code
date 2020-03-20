
function plotDistributionDensityOnlyOne(min_val,max_val,vector_one,label_one,title_label,x_label,y_label)

f = figure;

nPoints = 1000;

vector_one(vector_one==0) = [];
XI= linspace(min_val,max_val,nPoints);

kvone = ksdensity(vector_one,XI);
% kvone = kvone .* (sum(vector_one(:)));

idx = find(kvone==0);
kvone(idx) = [];
XIone = XI;
XIone(idx) = [];
plot(XIone,kvone,'r*');

xlabel(x_label);
ylabel(y_label);

min_y = min([kvone(:)]);
max_y = max([kvone(:)]);

box off
axis tight
legend(label_one);
legend boxoff

title(title_label);

xlim([(min_val) (max_val)]);
ylim([min_y  (max_y)]);

print(f,'-depsc',strcat('Distribution-',strrep(title_label,' ','-'),'.eps'));

end