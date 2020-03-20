
function plotDistributionDensity(vector_one,vector_two,label_one,label_two,title_label,x_label,y_label,color_one,color_two,plot_only_one,Should_I_Flip)


% vector_one = vector_one * scale;
% vector_two = vector_two * scale;

f = figure;

nPoints =  1000;

LineWidth = 7;

vector_one(vector_one==0) = [];
vector_two(vector_two==0) = [];
vector_one(isnan(vector_one)) = [];
vector_two(isnan(vector_two)) = [];
vector_one(isinf(vector_one)) = [];
vector_two(isinf(vector_two)) = [];

min_val = min([vector_one(:);vector_two(:)]);
max_val = max([vector_one(:);vector_two(:)]);

% min_val = 0;
% max_val = 1 * scale * 10;

XI = linspace(min_val,max_val,nPoints);

kvone = ksdensity(vector_one,XI);
kvtwo = ksdensity(vector_two,XI);
%kvone = kvone .* (sum(vector_one(:)));
%kvtwo = kvtwo .* (sum(vector_two(:)));

if Should_I_Flip; kvone = flip(kvone); kvtwo = flip(kvtwo); end

idx = find(kvone==0);
kvone(idx) = [];
XIone = XI;
XIone(idx) = [];
plot(XIone,kvone,color_one,'MarkerSize',1,'LineWidth',LineWidth);

hold on

if ~plot_only_one
    
    idx = find(kvtwo==0);
    kvtwo(idx) = [];
    XItwo = XI;
    XItwo(idx) = [];
    plot(XItwo,kvtwo,color_two,'MarkerSize',3,'LineWidth',LineWidth);

end

% xlabel(x_label);
% ylabel(y_label);

min_y = min([kvone(:); kvtwo(:)]);
max_y = max([kvone(:); kvtwo(:)]);

box off
axis tight
% if ~plot_only_one; legend(label_one,label_two); end
% legend boxoff

% title(title_label);

xlim([(min_val - 0.1) (max_val + 0.1)]);
% xlim([1 758]);
ylim([min_y max_y]);

set(gca,'YTick',[],'XTick',[]);

axis off

print(f,'-depsc',strcat('Distribution-',strrep(title_label,' ','-'),'.eps'));

end