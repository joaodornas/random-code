

function real_granger

% settings_subj1_2210;
% settings_subj2_2610;
% settings_subj3_0311;
% settings_subj4_0211;
% settings_subj5_0211;
% settings_subj6_2411;
% settings_subj7_1401;
% settings_subj8_1401;

%all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','DMN','AUD'};

%doEverythingAllSubjects;

%settings_subj8_1401_MAC;
%settings_subj7_1401_MAC;
%all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','AUD'};
%all_networks_labels = {'DAN','VAN','VIS'};
%doEverythingSubjectBySubjectOnMAC(settings);
%computeGrangerBtwNetsOnlySeedsMAC(all_networks_labels);

% settings_subj1_2210;
% settings_subj2_2610;
% settings_subj3_0311;
% settings_subj4_0211;
% settings_subj5_0211;
% settings_subj6_2411;
% all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','AUD'};
% all_networks_labels = {'DAN','VAN','VIS'};
% all_networks_labels = {'VAN','VIS'};
% computeGrangerBtwNetsOnlySeeds(settings,all_networks_labels);

settings_subj2_2610;
all_networks_labels = {'DAN','VAN','VIS'};
countSignificantGrangerPvals(settings,all_networks_labels);

end

function doEverythingAllSubjects

%only_functional_parcellation = {'DAN','VAN','SMN','VIS','FPC','LAN','DMN'};
% reduced_parcellation = {'DAN','VAN'};
%reduced_parcellation = {'DAN','VIS'};
all_networks = {'DAN','VIS','VAN','SMN','FPC','LAN','DMN','AUD'};

%nNet = length(reduced_parcellation);
nNet = length(all_networks);

all_settings = getAllSettings;

nSet = length(all_settings);

for iSet=5:6
    
    settings = all_settings(iSet).settings;
    
    disp(strcat(settings.codes.subject,'-',datestr(now)));

    computeGrangerBtwNetsOnlySeeds(settings,all_networks);
    
end

end

function doEverythingSubjectBySubjectOnMAC

%only_functional_parcellation = {'DAN','VAN','SMN','VIS','FPC','LAN','DMN'};
% reduced_parcellation = {'DAN','VAN'};
%reduced_parcellation = {'DAN','VIS'};
%all_networks = {'DAN','VIS','VAN','SMN','FPC','LAN','DMN','AUD'};
%all_networks = {'VIS','DAN','VAN','SMN','FPC','LAN','DMN','AUD'};
all_networks = {'DAN','VAN','SMN','VIS','FPC','AUD'};

%nNet = length(reduced_parcellation);
nNet = length(all_networks);

all_settings = getAllSettings;

nSet = length(all_settings);

for iSet=1:nSet
    
    settings = all_settings(iSet).settings;
    
    disp(strcat(settings.codes.subject,'-',datestr(now)));

    computeGrangerBtwNetsOnlySeedsMAC(settings,all_networks);
    
end

end

function computeGrangerBtwNets(settings,nets)

%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

real_load_all_data_FSL;

nNet = length(nets);

filename = strcat('LHR-',settings.codes.subject,'-Functional-Parcels-Pop-Map.nii');
foldername = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\','output','\','Functional_Parcellation','\','after_correction','\','Final_Parcellation');
load_parcellation = nifti(filename);
load_parcellation.dat.fname = strcat(foldername,'\',load_parcellation.dat.fname);

Functional_Parcellation_img = load_parcellation.dat(:,:,:);

analysis_label = 'Granger-Functional-Parcels';

for iNet=1:nNet
    
    network_label = nets{iNet};
    seeds = getFunctionalSeeds(network_label);
    
    all_seeds(iNet).network_label = network_label;
    
    nSeeds = length(seeds.ROI);
    
    for iSeed=1:nSeeds
        
       label_parcel = seeds.ROI(iSeed).label; 
       idx_parcel = seeds.ROI(iSeed).idx;
       
       all_seeds(iNet).seed(iSeed).label = label_parcel; 
       all_seeds(iNet).seed(iSeed).idx = idx_parcel; 
       
       idx_voxels = find(Functional_Parcellation_img==idx_parcel);
       
       all_seeds(iNet).seed(iSeed).idx_voxels = idx_voxels;
       
       all_seeds(iNet).seed(iSeed).nVoxels = length(idx_voxels);
        
    end 
    
end

for iNet=1:nNet
    
    nSeeds = length(all_seeds(iNet).seed);
    
    all_voxels_inet = [];
    
    for iSeed=1:nSeeds
        
        all_voxels_inet = [all_voxels_inet(:);all_seeds(iNet).seed(iSeed).idx_voxels(:)];
        
    end
    
    for iiNet=iNet:nNet

        disp(strcat('iNet:',all_seeds(iNet).network_label,'/','iiNet:',all_seeds(iiNet).network_label,'-',datestr(now)));
        
        nSeeds = length(all_seeds(iiNet).seed);

        all_voxels_iinet = [];

        for iSeed=1:nSeeds

            all_voxels_iinet = [all_voxels_iinet(:);all_seeds(iiNet).seed(iSeed).idx_voxels(:)];

        end    
        
        disp(strcat('iNetVoxels:',int2str(length(all_voxels_inet)),'/','iiNetVoxels:',int2str(length(all_voxels_iinet))));
        
        disp(strcat('Track','-',datestr(now)));
        [Nets_GF, Nets_Gpval] = computeGrangerAcrossVoxelsPerCondition(all_voxels_inet,all_voxels_iinet,Functional_Parcellation_img,Track);
       
        save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-','Track','-',all_seeds(iNet).network_label,'-',all_seeds(iiNet).network_label,'.mat'),'Nets_GF','Nets_Gpval','all_voxels_inet','all_voxels_iinet','all_seeds');
        
        disp(strcat('Passive','-',datestr(now)));
        [Nets_GF, Nets_Gpval] = computeGrangerAcrossVoxelsPerCondition(all_voxels_inet,all_voxels_iinet,Functional_Parcellation_img,Passive);
       
        save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-','Passive','-',all_seeds(iNet).network_label,'-',all_seeds(iiNet).network_label,'.mat'),'Nets_GF','Nets_Gpval','all_voxels_inet','all_voxels_iinet','all_seeds');
        
        disp(strcat('RestingState','-',datestr(now)));
        [Nets_GF, Nets_Gpval] = computeGrangerAcrossVoxelsPerCondition(all_voxels_inet,all_voxels_iinet,Functional_Parcellation_img,RestingState);
        
        save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-','RestingState','-',all_seeds(iNet).network_label,'-',all_seeds(iiNet).network_label,'.mat'),'Nets_GF','Nets_Gpval','all_voxels_inet','all_voxels_iinet','all_seeds');
        
    end
     
end

end

function [Nets_GF, Nets_Gpval] = computeGrangerAcrossVoxelsPerCondition(all_voxels_inet,all_voxels_iinet,Condition)
    
standard_img = zeros(91,109,91);

nRun = 4;

niNetVoxels = length(all_voxels_inet);
niiNetVoxels = length(all_voxels_iinet);

% if isequal(niNetVoxels,niiNetVoxels)
%     all_voxels = all_voxels_inet; 
% else
%     all_voxels = [all_voxels_inet;all_voxels_iinet]; 
% end

all_voxels = [all_voxels_inet;all_voxels_iinet];

nVoxels = length(all_voxels);

for iRun=1:nRun
    
    disp(strcat('Run:',int2str(iRun),'-',datestr(now)));

    for iVoxel=1:nVoxels
        
        if mod(iVoxel,1000) == 0, disp(int2str(iVoxel)); end

        [idxx,idxy,idxz] = ind2sub(size(standard_img),all_voxels(iVoxel));

        X = squeeze(Condition(iRun).run(idxx,idxy,idxz,:));
        
        for iiVoxel=iVoxel:nVoxels

            [iidxx,iidxy,iidxz] = ind2sub(size(standard_img),all_voxels(iiVoxel));
            
            Y = squeeze(Condition(iRun).run(iidxx,iidxy,iidxz,:));
            
            [GF,Gpval, GSig, morder, A, G, info] = getGrangerBtw2Voxels(X,Y);
            
            Nets_GF(iRun,iiVoxel,iVoxel) = GF(2,1);
            Nets_GF(iRun,iVoxel,iiVoxel) = GF(1,2);
            Nets_Gpval(iRun,iiVoxel,iVoxel) = Gpval(2,1);
            Nets_Gpval(iRun,iVoxel,iiVoxel) = Gpval(1,2);

        end

    end

end

end

function computeGrangerBtwNetsOnlySeeds(settings,nets)


%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

%real_load_all_data_FSL_without_trials;

nNet = length(nets);

analysis_label = 'Granger-Functional-Seeds';

for iNet=1:nNet
    
    network_label = nets{iNet};
    seeds = getFunctionalSeeds_v5(network_label);
    
    all_seeds(iNet).network_label = network_label;
    
    nSeeds = length(seeds.ROI);
    
    for iSeed=1:nSeeds
        
       label_parcel = seeds.ROI(iSeed).label; 
       idx_parcel = seeds.ROI(iSeed).idx;
       
       all_seeds(iNet).seed(iSeed).label = label_parcel; 
       all_seeds(iNet).seed(iSeed).idx = idx_parcel; 
       
       idx_voxels = getROISeedVoxelsIDX(seeds.ROI(iSeed));
       
       all_seeds(iNet).seed(iSeed).idx_voxels = idx_voxels;
       
       all_seeds(iNet).seed(iSeed).nVoxels = length(idx_voxels);
        
    end 
    
end

for iNet=1:nNet
    
    nSeeds = length(all_seeds(iNet).seed);
    
    all_voxels_inet = [];
    
    for iSeed=1:nSeeds
        
        all_voxels_inet = [all_voxels_inet(:);all_seeds(iNet).seed(iSeed).idx_voxels(:)];
        
    end
    
    for iiNet=iNet:nNet

        disp(strcat('iNet:',all_seeds(iNet).network_label,'/','iiNet:',all_seeds(iiNet).network_label,'-',datestr(now)));
        
        nSeeds = length(all_seeds(iiNet).seed);

        all_voxels_iinet = [];

        for iSeed=1:nSeeds

            all_voxels_iinet = [all_voxels_iinet(:);all_seeds(iiNet).seed(iSeed).idx_voxels(:)];

        end    
        
        disp(strcat('iNetVoxels:',int2str(length(all_voxels_inet)),'/','iiNetVoxels:',int2str(length(all_voxels_iinet))));
        
        kind = 'Track';
        for irun=1:4;
            [Track(irun).run, Track(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
        end

        disp(strcat('Track','-',datestr(now)));
        [Nets_GF, Nets_Gpval] = computeGrangerAcrossVoxelsPerCondition(all_voxels_inet,all_voxels_iinet,Track);
       
        save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-','Track','-',all_seeds(iNet).network_label,'-',all_seeds(iiNet).network_label,'.mat'),'Nets_GF','Nets_Gpval','all_voxels_inet','all_voxels_iinet','all_seeds');
        
        clear Track
        
        kind = 'Passive';
        for irun=1:4;
            [Passive(irun).run, Passive(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
        end

        disp(strcat('Passive','-',datestr(now)));
        [Nets_GF, Nets_Gpval] = computeGrangerAcrossVoxelsPerCondition(all_voxels_inet,all_voxels_iinet,Passive);
       
        save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-','Passive','-',all_seeds(iNet).network_label,'-',all_seeds(iiNet).network_label,'.mat'),'Nets_GF','Nets_Gpval','all_voxels_inet','all_voxels_iinet','all_seeds');
        
        clear Passive
        
        kind = 'RestingState';
        for irun=1:4;
            [RestingState(irun).run, RestingState(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
        end

        disp(strcat('RestingState','-',datestr(now)));
        [Nets_GF, Nets_Gpval] = computeGrangerAcrossVoxelsPerCondition(all_voxels_inet,all_voxels_iinet,RestingState);
        
        save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-','RestingState','-',all_seeds(iNet).network_label,'-',all_seeds(iiNet).network_label,'.mat'),'Nets_GF','Nets_Gpval','all_voxels_inet','all_voxels_iinet','all_seeds');
        
        clear RestingState
        
    end
     
end

end   

function idx_voxels = getROISeedVoxelsIDX(ROI)
        
MNI_size = [91, 109, 91];

MNI_x_center = 45;
MNI_y_center = 63;
MNI_z_center = 36;

size_voxels_mm = 2;
ROI_radius_mm = 6 - size_voxels_mm;
ROI_radius_voxels = ROI_radius_mm/size_voxels_mm;

x_center = MNI_x_center + round(ROI.x/size_voxels_mm) + 1;
y_center = MNI_y_center + round(ROI.y/size_voxels_mm) + 1;
z_center = MNI_z_center + round(ROI.z/size_voxels_mm) + 1;

xgv = (x_center-ROI_radius_voxels):(x_center+ROI_radius_voxels);
ygv = (y_center-ROI_radius_voxels):(y_center+ROI_radius_voxels);
zgv = (z_center-ROI_radius_voxels):(z_center+ROI_radius_voxels);

[X,Y,Z] = meshgrid(xgv,ygv,zgv);

idx = find(X);
nLine = size(X,1);

iVoxel = 0;

for iLine=1:nLine
    
    for iiLine=1:nLine
        
        for iiiLine=1:nLine
            
            iVoxel = iVoxel + 1;
   
            idx_voxels(iVoxel) = sub2ind(MNI_size,xgv(iLine),ygv(iiLine),zgv(iiiLine));
            
        end
        
    end
    
end

end

function computeGrangerBtwNetsOnlySeedsMAC(nets)

settings_subj8_1401_MAC;
%settings_subj7_1401_MAC;

%% LOAD DATA

get_at_this_preprocessed_step = settings.FSL.folders.custom;
file = settings.FSL.files.functional.custom.residual_voxel;
mask = settings.FSL.files.mask.custom;

%real_load_all_data_FSL_without_trials;

nNet = length(nets);

analysis_label = 'Granger-Functional-Seeds';

for iNet=1:nNet
    
    network_label = nets{iNet};
    seeds = getFunctionalSeeds_v5(network_label);
    
    all_seeds(iNet).network_label = network_label;
    
    nSeeds = length(seeds.ROI);
    
    for iSeed=1:nSeeds
        
       label_parcel = seeds.ROI(iSeed).label; 
       idx_parcel = seeds.ROI(iSeed).idx;
       
       all_seeds(iNet).seed(iSeed).label = label_parcel; 
       all_seeds(iNet).seed(iSeed).idx = idx_parcel; 
       
       idx_voxels = getROISeedVoxelsIDX(seeds.ROI(iSeed));
       
       all_seeds(iNet).seed(iSeed).idx_voxels = idx_voxels;
       
       all_seeds(iNet).seed(iSeed).nVoxels = length(idx_voxels);
        
    end 
    
end

for iNet=1:nNet
    
    nSeeds = length(all_seeds(iNet).seed);
    
    all_voxels_inet = [];
    
    for iSeed=1:nSeeds
        
        all_voxels_inet = [all_voxels_inet(:);all_seeds(iNet).seed(iSeed).idx_voxels(:)];
        
    end
    
    for iiNet=iNet:nNet

        disp(strcat('iNet:',all_seeds(iNet).network_label,'/','iiNet:',all_seeds(iiNet).network_label,'-',datestr(now)));
        
        nSeeds = length(all_seeds(iiNet).seed);

        all_voxels_iinet = [];

        for iSeed=1:nSeeds

            all_voxels_iinet = [all_voxels_iinet(:);all_seeds(iiNet).seed(iSeed).idx_voxels(:)];

        end    
        
        disp(strcat('iNetVoxels:',int2str(length(all_voxels_inet)),'/','iiNetVoxels:',int2str(length(all_voxels_iinet))));
        
        kind = 'Track';
        for irun=1:4;
            [Track(irun).run, Track(irun).mask, settings] = real_get_data_FSL_MAC(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
        end

        disp(strcat('Track','-',datestr(now)));
        [Nets_GF, Nets_Gpval] = computeGrangerAcrossVoxelsPerCondition(all_voxels_inet,all_voxels_iinet,Track);
       
        save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-','Track','-',all_seeds(iNet).network_label,'-',all_seeds(iiNet).network_label,'.mat'),'Nets_GF','Nets_Gpval','all_voxels_inet','all_voxels_iinet','all_seeds');
        
        clear Track
        
        kind = 'Passive';
        for irun=1:4;
            [Passive(irun).run, Passive(irun).mask, settings] = real_get_data_FSL_MAC(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
        end

        disp(strcat('Passive','-',datestr(now)));
        [Nets_GF, Nets_Gpval] = computeGrangerAcrossVoxelsPerCondition(all_voxels_inet,all_voxels_iinet,Passive);
       
        save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-','Passive','-',all_seeds(iNet).network_label,'-',all_seeds(iiNet).network_label,'.mat'),'Nets_GF','Nets_Gpval','all_voxels_inet','all_voxels_iinet','all_seeds');
        
        clear Passive
        
        kind = 'RestingState';
        for irun=1:4;
            [RestingState(irun).run, RestingState(irun).mask, settings] = real_get_data_FSL_MAC(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
        end

        disp(strcat('RestingState','-',datestr(now)));
        [Nets_GF, Nets_Gpval] = computeGrangerAcrossVoxelsPerCondition(all_voxels_inet,all_voxels_iinet,RestingState);
        
        save(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-','RestingState','-',all_seeds(iNet).network_label,'-',all_seeds(iiNet).network_label,'.mat'),'Nets_GF','Nets_Gpval','all_voxels_inet','all_voxels_iinet','all_seeds');
        
        clear RestingState
        
    end
     
end

end   

function countSignificantGrangerPvals(settings,nets)

nNet = length(nets);

analysis_label = 'Granger-Functional-Seeds';

pcriterion = 0.01;

for iRun=1:4
    
for iNet=1:nNet
    
    inetwork_label = nets{iNet};
    
    for iiNet=iNet:nNet
        
        iinetwork_label = nets{iiNet};
    
        load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-','Track','-',inetwork_label,'-',iinetwork_label,'.mat'));
        
        ninet = length(all_voxels_inet);
        niinet = length(all_voxels_iinet);
        
        if iNet == iiNet
            
            percentage = length(find(Nets_Gpval(iRun,:,:) < pcriterion)) / (((ninet)^2));
            
            results_track{iNet+1,1} = inetwork_label;
            results_track{1,iNet+1} = inetwork_label;
            
            results_track{iNet+1,iNet+1} = percentage;
        
        else
            
            percentage = length(find(Nets_Gpval(iRun,1:ninet,(ninet+1):end) < pcriterion)) / (((ninet*niinet)));
        
            results_track{iNet+1,1} = inetwork_label;
            results_track{1,iiNet+1} = iinetwork_label;
            
            results_track{iNet+1,iiNet+1} = percentage;
           
            percentage = length(find(Nets_Gpval(iRun,(ninet+1):end,1:ninet) < pcriterion)) / (((ninet*niinet)));
        
            results_track{iiNet+1,1} = iinetwork_label;
            results_track{1,iNet+1} = inetwork_label;
            
            results_track{iiNet+1,iNet+1} = percentage;
           
        end
        
        clear all_voxels_inet
        clear all_voxels_iinet
        clear Nets_GF
        clear Nets_Gpval
        clear all_seeds
        
        load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-','Passive','-',inetwork_label,'-',iinetwork_label,'.mat'));
        
                if iNet == iiNet
            
            percentage = length(find(Nets_Gpval(iRun,:,:) < pcriterion)) / (((ninet)^2));
            
            results_passive{iNet+1,1} = inetwork_label;
            results_passive{1,iNet+1} = inetwork_label;
            
            results_passive{iNet+1,iNet+1} = percentage;
        
        else
            
            percentage = length(find(Nets_Gpval(iRun,1:ninet,(ninet+1):end) < pcriterion)) / ((ninet*niinet));
        
            results_passive{iNet+1,1} = inetwork_label;
            results_passive{1,iiNet+1} = iinetwork_label;
            
            results_passive{iNet+1,iiNet+1} = percentage;
           
            percentage = length(find(Nets_Gpval(iRun,(ninet+1):end,1:ninet) < pcriterion)) / ((ninet*niinet));
        
            results_passive{iiNet+1,1} = iinetwork_label;
            results_passive{1,iNet+1} = inetwork_label;
            
            results_passive{iiNet+1,iNet+1} = percentage;
           
        end
        
        clear all_voxels_inet
        clear all_voxels_iinet
        clear Nets_GF
        clear Nets_Gpval
        clear all_seeds
        
        load(strcat(settings.codes.experiment,'-',settings.codes.subject,'-',analysis_label,'-','RestingState','-',inetwork_label,'-',iinetwork_label,'.mat'));
        
                if iNet == iiNet
            
            percentage = length(find(Nets_Gpval(iRun,:,:) < pcriterion)) / (((ninet)^2));
            
            results_rest{iNet+1,1} = inetwork_label;
            results_rest{1,iNet+1} = inetwork_label;
            
            results_rest{iNet+1,iNet+1} = percentage;
        
        else
            
            percentage = length(find(Nets_Gpval(iRun,1:ninet,(ninet+1):end) < pcriterion)) / (((ninet*niinet)));
        
            results_rest{iNet+1,1} = inetwork_label;
            results_rest{1,iiNet+1} = iinetwork_label;
            
            results_rest{iNet+1,iiNet+1} = percentage;
           
            percentage = length(find(Nets_Gpval(iRun,(ninet+1):end,1:ninet) < pcriterion)) / (((ninet*niinet)));
        
            results_rest{iiNet+1,1} = iinetwork_label;
            results_rest{1,iNet+1} = inetwork_label;
            
            results_rest{iiNet+1,iNet+1} = percentage;
           
        end
        
        clear all_voxels_inet
        clear all_voxels_iinet
        clear Nets_GF
        clear Nets_Gpval
        clear all_seeds
        
    end
    
end

results_run(iRun).results_track = results_track;
results_run(iRun).results_passive = results_passive;
results_run(iRun).results_rest = results_rest;

clear results_track
clear results_passive
clear results_rest

end

save('LHR-SUBJ2-Results-Causality.mat','results_run');

end

