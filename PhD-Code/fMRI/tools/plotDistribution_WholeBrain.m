function plotDistribution_WholeBrain

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


MD_xi = [];
PETRA_xi = [];
HCP_xi = [];
MD_xi_shifted = [];
PETRA_xi_shifted = [];
HCP_xi_shifted = [];
for iROI=1:90
    
    ROI_Label = AAL_ROI(iROI).Nom_L;
    ROI_Label = strrep(ROI_Label,'_','-');
    
    load(strcat('Distribution','-',KIND,'-',ROI_Label,'.mat'));
    
    MD(iROI).kv = kvone;
    MD(iROI).xi = XIone;
    MD(iROI).xi_shifted = int16( ceil( XIone.*100 ) ); 
    
    PETRA(iROI).kv = kvtwo;
    PETRA(iROI).xi = XItwo;
    PETRA(iROI).xi_shifted = int16( ceil( XItwo.*100 ) );
    
    HCP(iROI).kv = kvthree;
    HCP(iROI).xi = XIthree;
    HCP(iROI).xi_shifted = int16( ceil( XIthree.*100 ) );
    
    MD_xi = [MD_xi,XIone];
    PETRA_xi = [PETRA_xi,XItwo];
    HCP_xi = [HCP_xi,XIthree];
    
    MD_xi_shifted = [MD_xi_shifted,MD(iROI).xi_shifted];
    PETRA_xi_shifted = [PETRA_xi_shifted,PETRA(iROI).xi_shifted];
    HCP_xi_shifted = [HCP_xi_shifted,HCP(iROI).xi_shifted];
    
    clear kvone kvtwo kvthree XIone XItwo XIthree
    
end

min_MD_xi_shifted = min(MD_xi_shifted(:));
max_MD_xi_shifted = max(MD_xi_shifted(:));

min_PETRA_xi_shifted = min(PETRA_xi_shifted(:));
max_PETRA_xi_shifted = max(PETRA_xi_shifted(:));

min_HCP_xi_shifted = min(HCP_xi_shifted(:));
max_HCP_xi_shifted = max(HCP_xi_shifted(:));

all_MD_xi_shifted = min_MD_xi_shifted:1:max_MD_xi_shifted;
all_PETRA_xi_shifted = min_PETRA_xi_shifted:1:max_PETRA_xi_shifted;
all_HCP_xi_shifted = min_HCP_xi_shifted:1:max_HCP_xi_shifted;

all_MD_kv = zeros(size(all_MD_xi_shifted));
all_PETRA_kv = zeros(size(all_PETRA_xi_shifted));
all_HCP_kv = zeros(size(all_HCP_xi_shifted));

for iROI=1:nROIs
         
   nKV = length(MD(iROI).kv);
   
   for iKV=1:nKV
       
      all_MD_kv(find(all_MD_xi_shifted==MD(iROI).xi_shifted(iKV))) = all_MD_kv(find(all_MD_xi_shifted==MD(iROI).xi_shifted(iKV))) + MD(iROI).kv(iKV);
       
   end
   
   nKV = length(PETRA(iROI).kv);
   
   for iKV=1:nKV
       
      all_PETRA_kv(find(all_PETRA_xi_shifted==PETRA(iROI).xi_shifted(iKV))) = all_PETRA_kv(find(all_PETRA_xi_shifted==PETRA(iROI).xi_shifted(iKV))) + PETRA(iROI).kv(iKV);
       
   end
   
   nKV = length(HCP(iROI).kv);
   
   for iKV=1:nKV
       
      all_HCP_kv(find(all_HCP_xi_shifted==HCP(iROI).xi_shifted(iKV))) = all_HCP_kv(find(all_HCP_xi_shifted==HCP(iROI).xi_shifted(iKV))) + HCP(iROI).kv(iKV);
       
   end
   
end

f = figure;

all_MD_kv = all_MD_kv ./ nROIs;
all_PETRA_kv = all_PETRA_kv ./ nROIs;
all_HCP_kv = all_HCP_kv ./ nROIs;

all_MD_xi = ( all_MD_xi_shifted ) ;
all_PETRA_xi = ( all_PETRA_xi_shifted ) ;
all_HCP_xi = ( all_HCP_xi_shifted ) ;

idx = find(all_MD_kv==0);
all_MD_kv(idx) = [];
all_MD_xi(idx) = [];
plot(all_MD_xi,all_MD_kv,'r-','MarkerSize',3);

hold on

idx = find(all_PETRA_kv==0);
all_PETRA_kv(idx) = [];
all_PETRA_xi(idx) = [];
plot(all_PETRA_xi,all_PETRA_kv,'b-','MarkerSize',1);

hold on

idx = find(all_HCP_kv==0);
all_HCP_kv(idx) = [];
all_HCP_xi(idx) = [];
plot(all_HCP_xi,all_HCP_kv,'k-','MarkerSize',1);

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

print(f,'-depsc',strcat('Distribution-Mean-ROIs-',KIND,'.eps'));

end

