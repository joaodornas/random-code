for i=[1:4,6:12,14:92]
    %eval(['sf_cell' num2str(i) '=[sf_fitted_curves_200.fitted_curves.cell' num2str(i) ' sf_fitted_curves_800.fitted_curves.cell' num2str(i) '];']);
    eval(['tf_cell' num2str(i) '=[tf_fitted_curves_200.fitted_curves.cell' num2str(i) ' tf_fitted_curves_800.fitted_curves.cell' num2str(i) '];']);
end

clear tf_fitted_curves_200 tf_fitted_curves_800

mean_cells=zeros(size(tf_cell1));

for i=[1:4,6:12,14:92]
    tf_i=eval(['tf_cell' num2str(i)]);
    mean_cells=mean_cells+tf_i;
end

mean_cells=mean_cells./90;
mesh(1:1000,log2(0.05:0.05:32),mean_cells(1:640,:)); figure(gcf)

