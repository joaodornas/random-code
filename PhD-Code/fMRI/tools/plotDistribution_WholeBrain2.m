
function plotDistribution_WholeBrain2

KIND = 'Runs';
doTheMath( KIND );

KIND = 'Subjects';
doTheMath( KIND );

KIND = 'Average';
doTheMath( KIND );

end

function doTheMath( KIND )

nROIs = 90;

%% LOAD AAL
load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('/Users/joaodornas/Dropbox (joaodornas)/_Research/_toolBox/aal_for_SPM8/',load_aal.dat.fname);
AAL_img = load_aal.dat(:,:,:);
load_roi = load('ROI_MNI_V4_List.mat');
AAL_ROI = load_roi.ROI;


all_MD_xi = [];
all_PETRA_xi = [];
all_HCP_xi = [];
all_MD_kv = [];
all_PETRA_kv = [];
all_HCP_kv = [];
MD_xi_shifted = [];
PETRA_xi_shifted = [];
HCP_xi_shifted = [];
for iROI=1:90
    
    ROI_Label = AAL_ROI(iROI).Nom_L;
    ROI_Label = strrep(ROI_Label,'_','-');
    
    load(strcat('Distribution','-',KIND,'-',ROI_Label,'.mat'));
    
    MD(iROI).kv = kvone;
    MD(iROI).xi = XIone;
    MD(iROI).xi_shifted = int16( XIone.*10000 ); 
    
    PETRA(iROI).kv = kvtwo;
    PETRA(iROI).xi = XItwo;
    PETRA(iROI).xi_shifted = int16( XItwo.*10000 );
    
    HCP(iROI).kv = kvthree;
    HCP(iROI).xi = XIthree;
    HCP(iROI).xi_shifted = int16( XIthree.*10000 );
    
    all_MD_kv = [all_MD_kv,kvone];
    all_PETRA_kv = [all_PETRA_kv,kvtwo];
    all_HCP_kv = [all_HCP_kv,kvthree];
    
    all_MD_xi = [all_MD_xi,XIone];
    all_PETRA_xi = [all_PETRA_xi,XItwo];
    all_HCP_xi = [all_HCP_xi,XIthree];
    
    MD_xi_shifted = [MD_xi_shifted,MD(iROI).xi_shifted];
    PETRA_xi_shifted = [PETRA_xi_shifted,PETRA(iROI).xi_shifted];
    HCP_xi_shifted = [HCP_xi_shifted,HCP(iROI).xi_shifted];
    
    clear kvone kvtwo kvthree XIone XItwo XIthree
    
end

f = figure;

for iROI=1:nROIs
    
    kv = MD(iROI).kv;
    idx = find(kv==0);
    xi = MD(iROI).xi;
    kv(idx) = [];
    xi(idx) = [];
    plot(xi,kv,'r-','MarkerSize',3);

    hold on
    
    kv = PETRA(iROI).kv;
    idx = find(kv==0);
    xi = PETRA(iROI).xi;
    kv(idx) = [];
    xi(idx) = [];
    plot(xi,kv,'b-','MarkerSize',1);

    hold on

    kv = HCP(iROI).kv;
    idx = find(kv==0);
    xi = HCP(iROI).xi;
    kv(idx) = [];
    xi(idx) = [];
    plot(xi,kv,'k-','MarkerSize',1);

    if iROI==1
    
        ylabel('Probability');
        xlabel('Correlations (Fisher-transformed)');

        min_y = min([all_MD_kv(:); all_PETRA_kv(:); all_HCP_kv(:)]);
        max_y = max([all_MD_kv(:); all_PETRA_kv(:); all_HCP_kv(:)]);

        min_x = min([all_MD_xi(:); all_PETRA_xi(:); all_HCP_xi(:)]);
        max_x = max([all_MD_xi(:); all_PETRA_xi(:); all_HCP_xi(:)]);

        box off
        axis tight
        legend('Magdeburg','Petra','HCP');
        legend boxoff

        title('Distribution of Correlations');

        xlim([min_x max_x]);
        ylim([-2 max_y]);
    
    end

end

print(f,'-depsc',strcat('Distribution-All-ROIs-',KIND,'.eps'));

end

