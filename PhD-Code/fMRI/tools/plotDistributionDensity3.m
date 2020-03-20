
function plotDistributionDensity3(min_val,max_val,vector_one,vector_two,vector_three,label_one,label_two,label_three,title_label,x_label,y_label)

f = figure;

nPoints = 1000;

vector_one(vector_one==0) = [];
vector_two(vector_two==0) = [];
vector_three(vector_three==0) = [];
vector_one(isnan(vector_one)) = [];
vector_two(isnan(vector_two)) = [];
vector_three(isnan(vector_three)) = [];
vector_one(isinf(vector_one)) = [];
vector_two(isinf(vector_two)) = [];
vector_three(isinf(vector_three)) = [];

min_val = min([vector_one(:);vector_two(:);vector_three(:)]);
max_val = max([vector_one(:);vector_two(:);vector_three(:)]);

XI= linspace(min_val,max_val,nPoints);

kvone = ksdensity(vector_one,XI);
kvtwo = ksdensity(vector_two,XI);
kvthree = ksdensity(vector_three,XI);

idx = find(kvone==0);
kvone(idx) = [];
XIone = XI;
XIone(idx) = [];
plot(XIone,kvone,'r-','MarkerSize',3);

hold on

idx = find(kvtwo==0);
kvtwo(idx) = [];
XItwo = XI;
XItwo(idx) = [];
plot(XItwo,kvtwo,'b-','MarkerSize',1);

hold on

idx = find(kvthree==0);
kvthree(idx) = [];
XIthree = XI;
XIthree(idx) = [];
plot(XIthree,kvthree,'k-','MarkerSize',1);

xlabel(x_label);
ylabel(y_label);

min_y = min([kvone(:); kvtwo(:); kvthree(:)]);
max_y = max([kvone(:); kvtwo(:); kvthree(:)]);

box off
axis tight
legend(label_one,label_two,label_three);
legend boxoff

title(title_label);

xlim([(min_val - 0.1) (max_val + 0.1)]);
ylim([min_y max_y]);

print(f,'-depsc',strcat('Distribution-',strrep(title_label,' ','-'),'.eps'));

end