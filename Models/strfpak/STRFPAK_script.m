%  Hello and welcome to STRFPAK_script.  This file is a collection of
%  commands which mirror the instructions performed by the gui version of
%  STRFPAK for the cell in this output directory.
%  This code is dynamically generated as you click through STRFPAK.  
%
%  As of version 4.1, this code looks for two .mat files:
%  STRFPAK_script_parameters.mat and STRFPAK_script_dataset.mat.  The first
%  contains all your options and settings; the second contains the file
%  names of your input datasets and everything specific to that one cell.  
%  We hope it will be easy to "hack" STRFPAK_script_dataset.mat so that 
%  you can do computations for other cells and datasets without having to 
%  run the gui version of STRFPAK each time.  Good luck!
%  Here's an example of the kind of code which makes this script portable:
%  These lines of code assume you have loaded "STRFPAK_script_dataset.mat"
%  and that your original template data set was in
%  /Applications/MATLAB7/work/Theunissen/half_hash/STRFPAK_4.1/DemoData/,
%  and the current data is in pwd.  The current data happens to have the
%  same types of stim and resp names (which makes this easier on the
%  coder), but you get the idea (I hope).
%  for jj = 1:length(rawDS); rawDS{jj}.respfiles = ...
%  strrep(rawDS{jj}.respfiles, '/Applications/MATLAB7/work/Theunissen/half_hash/STRFPAK_4.1/DemoData/',[pwd '/']); end
%  for jj = 1:length(rawDS); rawDS{jj}.stimfiles = ...
%  strrep(rawDS{jj}.stimfiles, '/Applications/MATLAB7/work/Theunissen/half_hash/STRFPAK_4.1/DemoData/',[pwd '/']); end
%  Now make a global variable "outputPath" in pwd, save outputPath and
%  rawDS in the file "STRFPAK_script_dataset.mat" and copy
%  "STRFPAK_script_parameters.mat" and STRFPAK_script.m to the current
%  directory, and run STRFPAK_script.  I know, it's hard the first time,
%  but you can even automate this to suit your database, and soon it will
%  be no problem at all.  Enjoy!
clear global
load STRFPAK_script_parameters.mat %contains all your tol values, smoothing windows, etc.
load STRFPAK_script_dataset.mat %everything specific to the cell in question: filenames of datasets and the output directory.
running_in_script_mode = 'Yes'; % For those pesky statements where there's soemthing in the code you just have to modify.  In the GUI version, this variable doesn't exist.

%%%  NEW in version 4.4:  If STRFPAK typically guesses your input filenames
%%%  corectly, you can try commenting out "load STRFPAK_script_dataset.mat"
%%%  (You still need the line "load STRFPAK_script_parameters.mat")
%%%  and un-commenting the following:

% global rawDS outputPath
% rawDS = guess_rawDS;
% outputPath = fullfile(pwd,'Output');

%%%  End of automatic dataset detection.  Now, the super-lazy way to do
%%%  STRFs is:
%%%      % Run STRFPAK on a demo data set to make sure everything works, 
%%%      % Add the path of the Output directory of this sample cell to your MATLAB path 
%%%      % Make the changes in STRFPAK_script code comments (so that datasets are automatically detected) as above
%%%      % Move to a new data directory
%%%      % type STRFPAK_script


rec_make_dir(outputPath);

%%%^^^ begin preprocess


global rawDS DS
numfiles = length(rawDS);
global NBAND num_trials
tempWait = waitbar(0,...
    sprintf('Preprocessing 2-D movie data using Fourier Power Model'));
hashes_of_stims = {};
global the_checksum
the_checksum = '';
for ii=1:numfiles
    waitbar(ii/numfiles, tempWait);
    [path,name,ext,ver] = fileparts(rawDS{ii}.stimfiles);
    switch ext
        case {'.wav'}
            errordlg('This option is only for movie. ','Data Type Error','modal');
            close(tempWait);
            return;
        case {'.mat'}  % 2-D raw data (movie)
            stimstim = load(rawDS{ii}.stimfiles);
            flds = fieldnames(stimstim);
            if (length(flds) == 1)
                stimstim = getfield(stimstim, char(flds{1}));
            end

        case {'.dat', '.txt'}
            stimstim = dlmread(rawDS{ii}.stimfiles);

    end
    if length(size(stimstim))>= 3
        sDim = '2-D';
    else
        sDim = '1-D';
        errordlg('This option is only for movie. ','Data Type Error','modal');
        close(tempWait);
        return;
    end

    % Reshaping the 2-D movie into one big vector
    [xsize, ysize, framesize] = size(stimstim);
    global NBAND
    global startidx_v stopidx_v bwindow_v btakepower_v bzdc_v btempnl_v

    % Take fourier power of movie
    stim=do_cached_calc('movpower',stimstim, startidx_v, stopidx_v, bwindow_v,...
        btakepower_v, bzdc_v, btempnl_v);
    if ~isempty(the_checksum)
        hashes_of_stims{ii} = the_checksum; % calculated in do_cached_calc; the_checksum is a global variable
    end

    [NBAND, DS{ii}.nlen] = size(stim);

    save(fullfile(outputPath,[name,'_Stim_',num2str(ii),'.mat']), 'stim');
    if isempty(the_checksum)
        hashes_of_stims{ii} = checksum_from_file(fullfile(outputPath,[name,'_Stim_',num2str(ii),'.mat']));
    end
    % Assign values to global variable DS
    DS{ii}.stimfiles = fullfile(outputPath,[name,'_Stim_',num2str(ii),'.mat']);
    %DS{ii}.nlen = framesize;

    % 3.2. Take care of Response file
    %    If you have multiple trial data, calculate psth first
    %    Then resample it using new amp_samp_rate

    %rawResp = load(rawDS{ii}.respfiles);
    % Modified by Junli, 2003 to read new spike arrivial time file
    %
    [rpath,rname,rext,rver] = fileparts(rawDS{ii}.respfiles);
    switch rext
        case {'.dat', '.txt', ''}
            newpsth = load(rawDS{ii}.respfiles);

        case {'.mat'}
            respMat = load(rawDS{ii}.respfiles);

            % Validate the MAT-file
            flds = fieldnames(respMat);

            % Check if response data is in spike arrival time
            % or already preprocessed.
            if (length(flds) == 1)
                rawResp = getfield(respMat, char(flds{1}));
                if iscell(rawResp)
                    spiketrain = zeros(rawDS{ii}.ntrials,rawDS{ii}.nlen);
                    for trial_ind =1:rawDS{ii}.ntrials

                        spiketrain(trial_ind, rawResp{trial_ind}) = ones(1,length(rawResp{trial_ind}));
                    end
                    newpsth = resample(spiketrain', ampsamprate, respsamprate);

                    newpsth = newpsth'; % make sure new response data is trials x T.
                    newpsth(find(newpsth < 0)) = 0;
                else
                    newpsth = rawResp;
                end
                [xs, ys] = size(newpsth);
                if xs > ys   % make sure response data is trials X TimeFrame
                    newpsth = newpsth';
                end
                rgoodidx = find(~isnan(newpsth));
                newpsth = newpsth(rgoodidx);

            end
        otherwise
            errordlg(['Only ASCII and MAT binary fileformats work',...
                ' in this version.'], 'File Type Error', 'modal')
            return;
    end

    % save to the file for each data pair
    save(fullfile(outputPath,[name,'_Spike_time_',num2str(ii),'.mat']), 'newpsth');

    % Assign values to global variable DS
    DS{ii}.respfiles = fullfile(outputPath,[name,'_Spike_time_',num2str(ii),'.mat']);
    DS{ii}.ntrials = size(newpsth, 1);

end
put_stim_checksums(hashes_of_stims,DS);

global originalDS
originalDS = DS;
%%%^^^ end preprocess


%%%^^^ begin select_validation
%%%^^^ end select_validation

%%%^^^ begin calculate
%%%^^^ end calculate

%%%^^^ begin validation
%%%^^^ end validation

%%%^^^ begin goodness_of_fit
%%%^^^ end goodness_of_fit

%%%^^^ begin prediction
%%%^^^ end prediction


