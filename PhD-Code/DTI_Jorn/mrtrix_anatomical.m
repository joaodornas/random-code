
subject_list = {'SUBJECT-1-22-10-2015', 'SUBJECT-2-26-10-2015', 'SUBJECT-3-3-11-2015', 'SUBJECT-4-2-11-2015', 'SUBJECT-5-2-11-2015', 'SUBJECT-6-24-11-2015', 'SUBJECT-7-14-01-2016', 'SUBJECT-8-14-01-2016'};

for iSubj=1:1

  fileroot = strcat('/Volumes/dropbox/_DATA/LOW-HIGH-ATTENTION/',subject_list{iSubj});
  
  MPRAGEDIR = strcat(fileroot,'/preprocessed/T1-Anatomical/Run-1-0/FSL/original.bet');
  
  FILES = dir(MPRAGEDIR);
  MPRAGEFILE = FILES(end-2).name;
 
  system(sprintf('/Applications/mrtrix3/scripts/5ttgen ''fsl'' %s %s_ACT.nii',MPRAGEFILE,MPRAGEFILE(1:end-7)));
    
end