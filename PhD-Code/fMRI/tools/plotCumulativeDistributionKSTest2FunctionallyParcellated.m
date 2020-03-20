
function plotCumulativeDistributionKSTest2FunctionallyParcellated(matrix_one,matrix_two,label_one,label_two,title_label,x_label,y_label)


all_nets = getVoxelsPerNet;
nNet = 8;

%%% LOAD AAL

disp('Load AAL');
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('Z:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;

idx_AAL = find(AAL_img);

f = figure;

for iNet=1:nNet

    idx_AAL_Functional = idx_AAL(ismember(all_nets(iNet).idx_voxels,idx_AAL));
    idx_InSequence_AAL = translate3DAALtoSequence(idx_AAL_Functional);
    
    vector_one = matrix_one(

    vector_one(vector_one==0) = [];
    vector_two(vector_two==0) = [];

    [kvone xone] = ecdf(vector_one);
    [kvtwo xtwo] = ecdf(vector_two);

    idx_zeros_one = find(xone==0);
    idx_zeros_two = find(xtwo==0);
    kvone(idx_zeros_one) = [];
    xone(idx_zeros_one) = [];
    kvtwo(idx_zeros_two) = [];
    xtwo(idx_zeros_two) = [];

    [h,p] = kstest2(kvone,kvtwo);

    plot(xone,kvone,'ro','MarkerSize',6);
    hold on
    plot(xtwo,kvtwo,'b*','MarkerSize',2);

    x_text = 10000;
    y_text = 0.5;

    text(x_text,y_text,strcat('Kolmogorov-Smirnov test significance: p=',num2str(p)));

    xlabel(x_label);
    ylabel(y_label);

    box off
    axis tight
    legend(label_one,label_two);
    legend boxoff

    title(title_label);

    ylim([-0.5 1.5]);

end

print(f,'-depsc',strcat('CumulativeDistribution-KSTest2-',strrep(title_label,' ','-'),'.eps'));



end