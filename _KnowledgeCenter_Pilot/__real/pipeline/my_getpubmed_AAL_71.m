
function my_getpubmed_AAL_71

query_id = 1;
start_year = 1900;
end_year = 2016;
minutes_to_pause = 1;

ROIs = {'Precentral_L' 'Precentral_R' 'Frontal_Sup_L' 'Frontal_Sup_R' 'Frontal_Sup_Orb_L' 'Frontal_Sup_Orb_R' 'Frontal_Mid_L' 'Frontal_Mid_R' 'Frontal_Mid_Orb_L' 'Frontal_Mid_Orb_R' 'Frontal_Inf_Oper_L' 'Frontal_Inf_Oper_R' 'Frontal_Inf_Tri_L' 'Frontal_Inf_Tri_R' 'Frontal_Inf_Orb_L' 'Frontal_Inf_Orb_R' 'Rolandic_Oper_L' 'Rolandic_Oper_R' 'Supp_Motor_Area_L' 'Supp_Motor_Area_R' 'Olfactory_L' 'Olfactory_R' 'Frontal_Sup_Medial_L' 'Frontal_Sup_Medial_R' 'Frontal_Med_Orb_L' 'Frontal_Med_Orb_R' 'Rectus_L' 'Rectus_R' 'Insula_L' 'Insula_R' 'Cingulum_Ant_L' 'Cingulum_Ant_R' 'Cingulum_Mid_L' 'Cingulum_Mid_R' 'Cingulum_Post_L' 'Cingulum_Post_R' 'Hippocampus_L' 'Hippocampus_R' 'ParaHippocampal_L' 'ParaHippocampal_R' 'Amygdala_L' 'Amygdala_R' 'Calcarine_L' 'Calcarine_R' 'Cuneus_L' 'Cuneus_R' 'Lingual_L' 'Lingual_R' 'Occipital_Sup_L' 'Occipital_Sup_R' 'Occipital_Mid_L' 'Occipital_Mid_R' 'Occipital_Inf_L' 'Occipital_Inf_R' 'Fusiform_L' 'Fusiform_R' 'Postcentral_L' 'Postcentral_R' 'Parietal_Sup_L' 'Parietal_Sup_R' 'Parietal_Inf_L' 'Parietal_Inf_R' 'SupraMarginal_L' 'SupraMarginal_R' 'Angular_L' 'Angular_R' 'Precuneus_L' 'Precuneus_R' 'Paracentral_Lobule_L' 'Paracentral_Lobule_R' 'Caudate_L' 'Caudate_R' 'Putamen_L' 'Putamen_R' 'Pallidum_L' 'Pallidum_R' 'Thalamus_L' 'Thalamus_R' 'Heschl_L' 'Heschl_R' 'Temporal_Sup_L' 'Temporal_Sup_R' 'Temporal_Pole_Sup_L' 'Temporal_Pole_Sup_R' 'Temporal_Mid_L' 'Temporal_Mid_R' 'Temporal_Pole_Mid_L' 'Temporal_Pole_Mid_R' 'Temporal_Inf_L' 'Temporal_Inf_R'};

ROI_start = 71;
nROIs = 19;

for iROI=ROI_start:(ROI_start+nROIs)
    
    label = strcat('AAL-',ROIs{iROI});
   
    my_getpubmed(label,query_id,start_year,end_year,minutes_to_pause);
    
end

end