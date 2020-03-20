
function lowhigh_reslice_nii_viewer_zscore

    settings_jan_0502;
    
    fromZScoreToViewer(settings);

return

function fromZScoreToViewer(settings)

    output_zscore = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',settings.folders.output,'\','zscore');

    for run=1:2

        zscore_file = strcat(output_zscore,'\','Low-High-Jan-MOT4-Run-',int2str(run),'-zscore-mean-t.nii');
        viewer_file = strcat(output_zscore,'\','Low-High-Jan-MOT4-Run-',int2str(run),'-zscore-mean-t-viewer.nii');

        reslice_nii(zscore_file,viewer_file);

    end

    for run=1:2

        zscore_file = strcat(output_zscore,'\','Low-High-Jan-MOT2-Run-',int2str(run),'-zscore-mean-t.nii');
        viewer_file = strcat(output_zscore,'\','Low-High-Jan-MOT2-Run-',int2str(run),'-zscore-mean-t-viewer.nii');

        reslice_nii(zscore_file,viewer_file);

    end

    for run=1:2

        zscore_file = strcat(output_zscore,'\','Low-High-Jan-RestingState-Run-',int2str(run),'-zscore-mean-t.nii');
        viewer_file = strcat(output_zscore,'\','Low-High-Jan-RestingState-Run-',int2str(run),'-zscore-mean-t-viewer.nii');

        reslice_nii(zscore_file,viewer_file);

    end

return

