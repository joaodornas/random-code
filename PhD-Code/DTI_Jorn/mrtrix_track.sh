#!/usr/bin/csh -f

  mkdir $DATADIR/Fibres
  mkdir $DATADIR/Fibres/VP
  mkdir $DATADIR/Fibres/VP/NAc_2_Target_ACT

  mkdir $DATADIR/Fibres/VP/NAc_2_Target_ACT/NAc_Re_Target_Re_ACT

  set NUMBER = 100000
  echo $DATADIR
  echo number: $NUMBER

  # forward direction
  
  # NAc_2_Target_ACT/NAc_Re_Target_Re_ACT
  echo
  echo "NAc_2_Target_ACT/NAc_Re_Target_Re_ACT"
  echo
  set OUTPUTDIR = $DATADIR/Fibres/VP/NAc_2_Target_ACT/NAc_Re_Target_Re_ACT
  # input ROI finden
  set INPUT_ROI = `echo $T1ROIDIR/bin_0.5_coreg2diff_*_T1_Template_08-mask_NAc_Re.nii`
  echo $INPUT_ROI
  set len_NAc_2_VP = 10
  set EXCLUDE_ROI1 = $T1ROIDIR/Excl_ROI_N_opticus_distal.mif
  set EXCLUDE_ROI2 = $T1ROIDIR/Excl_ROI_N_opticus_proximal.mif
  set INCLUDE_ROI  = `echo $T1ROIDIR/bin_0.5_coreg2diff_*_T1_Template_08-mask_VP_Re.nii`
  tckgen $DATADIR/SH_FOD_file.nii $OUTPUTDIR/NAc_Re_to_VP_Re_b1000_"$len_NAc_2_VP"_mm_"$SID"_ACT.tck -seed_image $INPUT_ROI -include $INCLUDE_ROI -act $ACT_FILE -exclude $EXCLUDE_ROI1 -exclude $EXCLUDE_ROI2 -force -number $NUMBER -fslgrad $DATADIR/bvecs $DATADIR/bvals -nthreads 4 -info -unidirectional -maxlength $len_NAc_2_VP

  #creation of density maps (the same as fdt_path in fsl)
  cd $OUTPUTDIR
  foreach map_tck (*ACT.tck)
    echo $map_tck
    tckmap $OUTPUTDIR/$map_tck -template $MPRAGE_2_Diff $OUTPUTDIR/"map_16_ACT_"$map_tck.nii -force
    tckmap $OUTPUTDIR/$map_tck -template {$MPRAGEFILE}.nii $OUTPUTDIR/"map_08_ACT_"$map_tck.nii -force
  end