

%file = settings.FSL.files.functional.warped;
%file = settings.FSL.files.functional.masked;
%file = settings.FSL.files.functional.despike;

%file = settings.FSL.files.functional.custom.residual_voxel;
%mask = settings.FSL.files.mask.custom;

%mask = settings.FSL.files.mask.warped;
%mask = settings.FSL.files.mask.custom;

%get_at_this_preprocessed_step = settings.FSL.folders.warped;
%get_at_this_preprocessed_step = settings.FSL.folders.custom;

disp(file);

kind = 'Track';
for irun=1:4;
    [Track(irun).run, Track(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
end

kind = 'Passive';
for irun=1:4;
    [Passive(irun).run, Passive(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
end

% kind = 'Trials';
% for irun=1:4;
%     [Trials(irun).run, Trials(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
% end

kind = 'RestingState';
for irun=1:4;
    [RestingState(irun).run, RestingState(irun).mask, settings] = real_get_data_FSL(settings,kind,irun,file,mask,get_at_this_preprocessed_step);
end
