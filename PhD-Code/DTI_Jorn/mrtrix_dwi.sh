#!/usr/bin/csh -f

  #create mask.nii
  /Applications/mrtrix3/release/bin/dwi2mask -info -fslgrad 6.forbedpost/bvecs 6.forbedpost/bvals -force -debug 6.forbedpost/data.nii mrtrix/mrtrix_mask

  #compute response funktion
  dwi2response msmt_5tt 6.forbedpost/data.nii.gz ../T1-Anatomical/Run-1-0/FSL/original.bet/co20151026_120654T1MPRAGE1x1x1sagp2s005a1001_ACT.nii.gz mrtrix/out_wm.txt mrtrix/out_gm.txt mrtrix/out_csf.txt -shell 0,2000 -lmax 10,10 -mask mrtrix/mrtrix_mask.nii -voxels mrtrix/single-fibre-voxels.nii -fslgrad 6.forbedpost/bvecs 6.forbedpost/bvals -fa 0.15

  #fit the SH-FOD-file (contrained spherical harmonics)
<<<<<<< Local Changes
  /Applications/mrtrix3/release/bin/dwi2fod -fslgrad 6.forbedpost/bvecs 6.forbedpost/bvals -shell 0,2000 -mask mrtrix3/mrtrix_mask.nii -lmax 10,10 -nthreads 4 -force -info 6.forbedpost/data.nii mrtrix/response_out mrtrix3/SH_FOD_file
=======
  /Applications/mrtrix3/release/bin/dwi2fod -fslgrad 6.forbedpost/bvecs 6.forbedpost/bvals -shell 0,2000 -mask mrtrix/mrtrix_mask.nii -lmax 10,10 -nthreads 4 -force -info msmt_csd 6.forbedpost/data.nii.gz mrtrix/out_wm.txt mrtrix/out_wm.mif mrtrix/out_gm.txt mrtrix/out_gm.mif mrtrix/out_csf.txt mrtrix/out_csf.mif
>>>>>>> External Changes



