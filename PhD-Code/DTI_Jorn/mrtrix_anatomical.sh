#!/usr/bin/csh -f

# subject list
#set subject_list = "SUBJECT-1-22-10-2015 SUBJECT-2-26-10-2015 SUBJECT-3-3-11-2015 SUBJECT-4-2-11-2015 SUBJECT-5-2-11-2015 SUBJECT-6-24-11-2015 SUBJECT-7-14-01-2016 SUBJECT-8-14-01-2016"      

set subject_list = "SUBJECT-1-22-10-2015"         
echo $subject_list

foreach SID ($subject_list)                      #if you have a subject list
 
  echo $SID
  echo
 
  set fileroot = /Volumes/dropbox/_DATA/LOW-HIGH-ATTENTION/$SID  #an example folder
  set MPRAGEDIR = $fileroot/preprocessed/T1-Anatomical/Run-1-0/FSL/original.bet   # an subfolder with MPRAGE data
  # find MPRAGE file
  set MPRAGEFILE = `echo {$MPRAGEDIR}/co*MPRAGE*001.nii.gz | awk '{print(substr($1,1, length($1)-4))}'`  # MPRAGE file is converted with dcm2niigui (mricron) from DICOM
 
  echo "MPRAGE:"
  echo $MPRAGEFILE
  echo " "
 
  # perform ACT
  /Applications/mrtrix3/scripts/5ttgen {$MPRAGEFILE}.nii.gz {$MPRAGEFILE}_ACT.nii -nthreads 4 fsl
    
end #over subject list
