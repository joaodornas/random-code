
function plotDistributionDensityXInLog(min_val,max_val,vector_one,vector_two,label_one,label_two,title_label,x_label,y_label)

f = figure;

nPoints = 1000;

vector_one(vector_one==0) = [];
vector_two(vector_two==0) = [];
XI= linspace(min_val,max_val,nPoints);

kvone = ksdensity(vector_one,XI);
kvtwo = ksdensity(vector_two,XI);
%kvone = kvone .* (sum(vector_one(:)));
%kvtwo = kvtwo .* (sum(vector_two(:)));

idx = find(kvone==0);
kvone(idx) = [];
XIone = XI;
XIone(idx) = [];
semilogx(XIone,kvone,'r-','MarkerSize',3);

hold on

idx = find(kvtwo==0);
kvtwo(idx) = [];
XItwo = XI;
XItwo(idx) = [];
semilogx(XItwo,kvtwo,'b-','MarkerSize',1);

xlabel(x_label);
ylabel(y_label);

min_y = min([kvone(:); kvtwo(:)]);
max_y = max([kvone(:); kvtwo(:)]);

box off
axis tight
legend(label_one,label_two);
legend boxoff

title(title_label);

xlim([(min_val) (max_val)]);
ylim([-0.0001  0.001]);

print(f,'-depsc',strcat('Distribution-',strrep(title_label,' ','-'),'.eps'));

end