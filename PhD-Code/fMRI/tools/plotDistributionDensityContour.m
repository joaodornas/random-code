function plotDistributionDensityContour(vector_one,vector_two,title_label,label_one,label_two)

f = figure;

vector_one(vector_one==0) = [];
vector_two(vector_two==0) = [];

Y_one = quantile(vector_one,[0.0025:0.0025:0.9975]);
Y_two = quantile(vector_two,[0.0025:0.0025:0.9975]);

[X,Y] = meshgrid(Y_one,Y_two);

Y_three = X ./ Y;

contour(Y_one,Y_two,Y_three,'ShowText','on');

box off
axis tight

xlim([min(Y_one) max(Y_one)]);
ylim([min(Y_two) max(Y_two)]);

xlabel(label_one);
ylabel(label_two);

title(title_label);

print(f,'-depsc',strcat('Contour-',strrep(title_label,' ','-'),'.eps'));

end