function lowhigh_Plot_3D_network_AAL_FC_DTI

settings_jan_0805;
%settings_elena_2905;

DTI_datafile = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',settings.folders.preprocessed,'\',settings.folders.DTI.main,'\',settings.folders.DTI.bedpostx,'\',settings.DTI.connectivity_matrix,'.mat');

% load your AAL network named C
con_mat = load(DTI_datafile);
C = con_mat.C;

load aal_cog.txt % Only 90 areas (not 116)
MNI_coord=aal_cog;
clear aal_cog

N=length(C);

figure('name','Connectome DIFF 3D')
set(gcf, 'color', 'white');

hold on

cmax=max(C(:));
C=(C/cmax); % In this way the strongest connection is =1.

% PLOT THE LINKS AS 3D CILINDERS
for n=1:N-1
    for p = n+1:N
        if C(n,p) > 0.80 % The cylinder's diameter is scaled proportionally to C
             [X,Y,Z]=cylinder1(MNI_coord(n,:),MNI_coord(p,:),C(n,p),10);
             surf(X,Y,Z,'FaceColor','b','EdgeColor','none') % Color in [Red Green Blue]
        elseif C(n,p) > 0.05 % between 0.15 and 0.01 all links have the minimun diameter of 0.6 
            [X,Y,Z]=cylinder1(MNI_coord(n,:),MNI_coord(p,:),0.2,10);
            surf(X,Y,Z,'FaceColor','b','EdgeColor','none') 
        end % The links below 0.01 (i.e less than 1% of the strongest link) are not plotted.
    end
end
% You can also scale the face color as a function of C(n,p).

% PLOT THE NODES AS 3D SPHERES
[x,y,z] = sphere;
for n=1:N
    surf(3*x+MNI_coord(n,1), 3*y+MNI_coord(n,2),3*z+MNI_coord(n,3),'FaceColor','r','EdgeColor','none');
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
D=spm_eeg_load;
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