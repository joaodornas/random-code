
function lowhigh_FC_AAL_DTI_3D

%%% SETTINGS

%settings_jan_0805;
settings_elena_2905;

%filename = 'filtered';
filename = 'residual';
%doTheMath(settings,filename);

plotResults(settings,filename);

end

function doTheMath(settings,filename)

%get_at_this_preprocessed_step = settings.FSL.folders.warped;
get_at_this_preprocessed_step = settings.FSL.folders.custom;
%file = settings.FSL.files.functional.custom.filtered;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

%% LOAD DATA

lowhigh_load_all_data_FSL;

%%% LOAD AAL

load_aal = nifti('ROI_MNI_V4.nii');
load_aal.dat.fname = strcat('K:\Dropbox (Uni Magdeburg)\_TOOLBOX\aal_for_SPM8\',load_aal.dat.fname);

AAL_img = load_aal.dat(:,:,:);

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

nNodes = 90;

%%% GET CONDITIONS

% MOT4 = cat(4,MOT4Run1,MOT4Run2);
% MOT2 = cat(4,MOT2Run1,MOT2Run2);
% RestingState = cat(4,RestingStateRun1,RestingStateRun2);

MOT4 = ( MOT4Run1 + MOT4Run2 ) ./2;
MOT2 = ( MOT2Run1 + MOT2Run2 ) ./2;
RestingState = ( RestingStateRun1 + RestingStateRun2 ) ./2;

% MOT4 = MOT4Run1;
% MOT2 = MOT2Run1;
% RestingState = RestingStateRun1;

%%% GET ROI MEAN TIME SERIES

nTR = size(MOT4,4);

mean_area_MOT4_voxel = zeros(nNodes,nTR);
mean_area_MOT2_voxel = zeros(nNodes,nTR);
mean_area_RestingState_voxel = zeros(nNodes,nTR);

%ROI_Zeros = zeros(nNodes,4);

for iNode=1:nNodes

    idx_ROI = AAL_ROI(iNode).ID;
    idx_ROI_Label = AAL_ROI(iNode).Nom_L;
    idx_ROI_Label = strrep(idx_ROI_Label,'_','-');
    
    idx_voxels = find(AAL_img == idx_ROI);
    
    nVoxels = length(idx_voxels);
    
    area_MOT4 = zeros(nVoxels,nTR);
    area_MOT2 = zeros(nVoxels,nTR);
    area_RestingState = zeros(nVoxels,nTR);
    
    for iVoxel=1:nVoxels
       
        [idxx,idxy,idxz] = ind2sub(size(AAL_img),idx_voxels(iVoxel));
            
        area_MOT4(iVoxel,:) = MOT4(idxx,idxy,idxz,:);
        
        area_MOT2(iVoxel,:) = MOT2(idxx,idxy,idxz,:);
            
        area_RestingState(iVoxel,:) = RestingState(idxx,idxy,idxz,:);

    end
    
    mn_MOT4 = mean(area_MOT4,1);
    mn_MOT2 = mean(area_MOT2,1);
    mn_RestingState = mean(area_RestingState,1);
      
    mean_area_MOT4_voxel(iNode,:) = mn_MOT4(:);
    mean_area_MOT2_voxel(iNode,:) = mn_MOT2(:);
    mean_area_RestingState_voxel(iNode,:) = mn_RestingState(:);
    
end


%%% GET CORRELATIONS

% mn_MOT4 = mean(mean_area_MOT4_voxel,2);                % compute row mean
% mean_area_MOT4_voxel = mean_area_MOT4_voxel - repmat(mn_MOT4,1,nTR);         % subtract row mean
% 
% mn_MOT2 = mean(mean_area_MOT2_voxel,2);                % compute row mean
% mean_area_MOT2_voxel = mean_area_MOT2_voxel - repmat(mn_MOT2,1,nTR);         % subtract row mean
% 
% mn_RestingState = mean(mean_area_RestingState_voxel,2);                % compute row mean
% mean_area_RestingState_voxel = mean_area_RestingState_voxel - repmat(mn_RestingState,1,nTR);         % subtract row mean

mean_area_MOT4_voxel = mean_area_MOT4_voxel';
mean_area_MOT2_voxel = mean_area_MOT2_voxel';
mean_area_RestingState_voxel = mean_area_RestingState_voxel';

[rho_MOT4,pval_MOT4] = corr(mean_area_MOT4_voxel);
[rho_MOT2,pval_MOT2] = corr(mean_area_MOT2_voxel);
[rho_RestingState,pval_RestingState] = corr(mean_area_RestingState_voxel);

increased_MOT4_MOT2 = rho_MOT4 > rho_MOT2;
increased_MOT2_RestingState = rho_MOT2 > rho_RestingState;

linear_increased_MOT4_MOT2_RestingState = (increased_MOT4_MOT2 == 1) & (increased_MOT2_RestingState);

k = 0;

for iNode=1:nNodes
    
    for iiNode=1:nNodes
        
        if linear_increased_MOT4_MOT2_RestingState(iNode,iiNode) == 1

            k = k + 1;
            
            linear_increased_labels{k,1} = AAL_ROI(iNode).Nom_L;
            linear_increased_labels{k,2} = AAL_ROI(iiNode).Nom_L;
            linear_increased_labels{k,3} = rho_MOT4(iNode,iiNode);
            linear_increased_labels{k,4} = rho_MOT2(iNode,iiNode);
            linear_increased_labels{k,5} = rho_RestingState(iNode,iiNode);

        end
        
    end
    
end

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-DTI-3D-',filename,'.mat'),'rho_MOT4','pval_MOT4','rho_MOT2','pval_MOT2','rho_RestingState','pval_RestingState','increased_MOT4_MOT2','increased_MOT2_RestingState','linear_increased_MOT4_MOT2_RestingState','linear_increased_labels','-v7.3');

end

function plotResults(settings,filename)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-DTI-3D-',filename,'.mat'));

load_roi = load('ROI_MNI_V4_List.mat');

AAL_ROI = load_roi.ROI;

%%% PLOT FC Matrix

CLOW = min([min(min(rho_MOT4)), min(min(rho_MOT2)), min(min(rho_RestingState))]);
CHIGH = max([max(max(rho_MOT4)), max(max(rho_MOT2)), max(max(rho_RestingState))]);

CLIM = [CLOW CHIGH];

f = figure;

imagesc(rho_MOT4,CLIM);
colorbar;
title('FC matrix');
xlabel('AAL regions');
ylabel('AAL regions');

print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-HighAttention','.jpg'));
print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-HighAttention','.eps'));
print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-HighAttention','.pdf'));

g = figure;

imagesc(rho_MOT2,CLIM);
colorbar;
title('FC matrix');
xlabel('AAL regions');
ylabel('AAL regions');

print(g,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-LowAttention','.jpg'));
print(g,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-LowAttention','.eps'));
print(g,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-LowAttention','.pdf'));

h = figure;

imagesc(rho_RestingState,CLIM);
colorbar;
title('FC matrix');
xlabel('AAL regions');
ylabel('AAL regions');

print(h,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-RestingState','.jpg'));
print(h,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-RestingState','.eps'));
print(h,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','FC-matrix-RestingState','.pdf'));

%%% PLOT 3D Rendering

DTI_datafile = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',settings.folders.preprocessed,'\',settings.folders.DTI.main,'\',settings.folders.DTI.bedpostx,'\',settings.DTI.connectivity_matrix,'.mat');

% load your AAL network named C
con_mat = load(DTI_datafile);
C = con_mat.C;

FC = linear_increased_MOT4_MOT2_RestingState;
FC = FC + 0;

load aal_cog.txt % Only 90 areas (not 116)
MNI_coord=aal_cog;
clear aal_cog

N=length(C);

figure('name','Connectome DIFF 3D')
set(gcf, 'color', 'white');

hold on

cmax=max(C(:));
C=(C/cmax); % In this way the strongest connection is =1.

Connections = applyFisherTransform(rho_MOT4,rho_MOT2,rho_RestingState,C);

% Connections = applyMcNemar(rho_MOT4,rho_MOT2,rho_RestingState,C);

% Connections = (C > 0.05) & (FC == 1) & (pval_MOT4 < 0.05);

Connections = Connections + 0;

FontSize = 12;

% PLOT THE LINKS AS 3D CILINDERS
for n=1:N-1
    for p = n+1:N
        if Connections(n,p) > 0.80 % The cylinder's diameter is scaled proportionally to C
             [X,Y,Z]=cylinder1(MNI_coord(n,:),MNI_coord(p,:),Connections(n,p),10);
             surf(X,Y,Z,'FaceColor','b','EdgeColor','none') % Color in [Red Green Blue]
        elseif Connections(n,p) > 0.05 % between 0.15 and 0.01 all links have the minimun diameter of 0.6 
            [X,Y,Z]=cylinder1(MNI_coord(n,:),MNI_coord(p,:),0.2,10);
            surf(X,Y,Z,'FaceColor','b','EdgeColor','none') 
        end % The links below 0.01 (i.e less than 1% of the strongest link) are not plotted.
    end
end
% You can also scale the face color as a function of C(n,p).

% PLOT THE NODES AS 3D SPHERES
nNode = size(Connections,1);
nodes = [];
for iNode=1:nNode

   if ~isempty(find(Connections(iNode,:)))
       
       nodes = [nodes, iNode];
       
   end

end
[x,y,z] = sphere;
%for n=1:N
for n=nodes
    surf(3*x+MNI_coord(n,1), 3*y+MNI_coord(n,2),3*z+MNI_coord(n,3),'FaceColor','r','EdgeColor','none');
    text(MNI_coord(n,1),MNI_coord(n,2),max(max(3*z+MNI_coord(n,3))),strrep(AAL_ROI(n).Nom_L,'_','-'),'FontSize',FontSize);
end
axis off;
axis equal

camlight;
rotate3d;
material dull; 
lighting phong;

%----------------- To plot the cortex/scalp underneath-------------------------

% put your own spm path 
%addpath spm8
%addpath spm8/canonical
D=spm_eeg_load('Neuromag306_EEG70.mat');
% Only forthis template
fn = fieldnames(D.inv{1}.mesh);
tess_ind = find(strncmp(fn,'tess_',5) & ~strcmp(fn,'tess_mni'));
for i = 1:length(tess_ind)
    location=D.inv{1}.mesh.(fn{tess_ind(i)});
    D.inv{1}.mesh.(fn{tess_ind(i)})=location(14:end);
end
mesh = spm_eeg_inv_transform_mesh(D.inv{1}.mesh.Affine, D.inv{1}.mesh);
Mscalp  = mesh.tess_scalp;
headface    = Mscalp.face;
headvert    = Mscalp.vert;
cortexface  = mesh.tess_ctx.face;
cortexvert  = mesh.tess_ctx.vert;


%plot3(headface(:,1),headface(:,2),headface(:,3),'.','MarkerSize',1,'Color',[0.95 .7 .55])
%patch('vertices',headvert,'faces',headface,'EdgeColor','none','FaceColor',[1 .7 .55],'FaceAlpha',0.1); %plots a mesh surface for the headshape
patch('vertices',cortexvert,'faces',cortexface,'EdgeColor','none','FaceColor',[0.5 0.5 .5],'FaceAlpha',0.15); %plots a mesh surface for the cortex

camlight;
rotate3d;
material dull; 
lighting phong;
%}

end

function Connections = applyFisherTransform(rho_MOT4,rho_MOT2,rho_RestingState,C)

fisher_MOT4 = (0.5) .* log( ( 1 + rho_MOT4 ) ./ ( 1 - rho_MOT4) );
fisher_MOT2 = (0.5) .* log( ( 1 + rho_MOT2 ) ./ ( 1 - rho_MOT2) );
fisher_RestingState = (0.5) .* log( ( 1 + rho_RestingState ) ./ ( 1 - rho_RestingState) );

nTR = 331;
n = nTR;

z_test_MOT4_MOT2 = ( fisher_MOT4 - fisher_MOT2 ) ./ sqrt( ( 1/(n-3) ) + ( 1/(n-3) ) );
z_test_MOT2_RestingState = ( fisher_MOT2 - fisher_RestingState ) ./ sqrt( ( 1/(n-3) ) + ( 1/(n-3) ) );

onetailed_MOT4_MOT2 = 1-normcdf(z_test_MOT4_MOT2,0,1);
onetailed_MOT2_RestingState = 1-normcdf(z_test_MOT2_RestingState,0,1);

Connections = (C > 0.05) & (onetailed_MOT4_MOT2 < 0.05) & (onetailed_MOT2_RestingState < 0.05);

end

function Connections = applyMcNemar(rho_MOT4,rho_MOT2,rho_RestingState,C)

chi_MOT4_MOT2 = ( rho_MOT4 - rho_MOT2 ).^2 ./ ( rho_MOT4 + rho_MOT2 ) ;

chi_MOT2_RestingState = ( rho_MOT2 - rho_RestingState ).^2 ./ ( rho_MOT2 + rho_RestingState ) ;

pval_MOT4_MOT2 = 1-chi2cdf(chi_MOT4_MOT2,1);

pval_MOT2_RestingState = 1-chi2cdf(chi_MOT2_RestingState,1);

Connections = (C > 0.05) & (pval_MOT4_MOT2 < 0.05) & (pval_MOT2_RestingState < 0.05);

end