
function roi_files_cortex_RL = get_ROI_aal_files_cortex_RL

toolbox_folder = 'K:\Dropbox (Uni Magdeburg)\_TOOLBOX\';
marsbar_aal_folder = 'marsbar-0.44\marsbar-aal-0.2\';
sufix = '_roi.mat';

ROI = getROI;

roi_files_cortex_RL = char.empty;
for i=1:length(ROI)
    
    roi_string = [];
    roi_string = strcat(toolbox_folder,marsbar_aal_folder,ROI{i},sufix);
    
    roi_files_cortex_RL = strvcat(roi_files_cortex_RL,roi_string);
    
end

end

function ROI = getROI
                           
    ROI{1} = 'MNI_Precentral_R';  
    ROI{2} = 'MNI_Paracentral_Lobule_R';  
    ROI{3} = 'MNI_Postcentral_R'; 
    ROI{4} = 'MNI_Frontal_Inf_Oper_R';  
    ROI{5} = 'MNI_Frontal_Inf_Orb_R';   
    ROI{6} = 'MNI_Frontal_Inf_Tri_R';   
    ROI{7} = 'MNI_Frontal_Med_Orb_R';   
    ROI{8} = 'MNI_Frontal_Mid_Orb_R';   
    ROI{9} = 'MNI_Frontal_Mid_R';       
    ROI{10} = 'MNI_Frontal_Sup_Medial_R';
    ROI{11} = 'MNI_Frontal_Sup_Orb_R';   
    ROI{12} = 'MNI_Frontal_Sup_R';     
    ROI{13} = 'MNI_Rolandic_Oper_R'; 
    ROI{14} = 'MNI_Supp_Motor_Area_R';
    
    ROI{15} = 'MNI_Parietal_Inf_R';         
    ROI{16} = 'MNI_Parietal_Sup_R'; 
    ROI{17} = 'MNI_Rectus_R';
    ROI{18} = 'MNI_ParaHippocampal_R';  
    ROI{19} = 'MNI_SupraMarginal_R'; 

    ROI{20} = 'MNI_Precuneus_R';
    ROI{21} = 'MNI_Angular_R';
    ROI{22} = 'MNI_Cuneus_R';            
    ROI{23} = 'MNI_Calcarine_R';    
    ROI{24} = 'MNI_Occipital_Inf_R';      
    ROI{25} = 'MNI_Occipital_Mid_R';        
    ROI{26} = 'MNI_Occipital_Sup_R';       
    
    ROI{27} = 'MNI_Olfactory_R'; 
    ROI{28} = 'MNI_Lingual_R';     
    ROI{29} = 'MNI_Heschl_R';
    ROI{30} = 'MNI_Fusiform_R';                
    ROI{31} = 'MNI_Insula_R';   
    ROI{32} = 'MNI_Temporal_Inf_R';          
    ROI{33} = 'MNI_Temporal_Mid_R';      
    ROI{34} = 'MNI_Temporal_Pole_Mid_R'; 
    ROI{35} = 'MNI_Temporal_Pole_Sup_R';    
    ROI{36} = 'MNI_Temporal_Sup_R';                  
    ROI{37} = 'MNI_Putamen_R';                    
      
%     ROI{2} = 'MNI_Amygdala_R';     
%     ROI{8} = 'MNI_Caudate_R';           
%     ROI{10} = 'MNI_Cerebelum_10_R';       
%     ROI{12} = 'MNI_Cerebelum_3_R';       
%     ROI{14} = 'MNI_Cerebelum_4_5_R';     
%     ROI{16} = 'MNI_Cerebelum_6_R';       
%     ROI{18} = 'MNI_Cerebelum_7b_R';      
%     ROI{20} = 'MNI_Cerebelum_8_R';       
%     ROI{22} = 'MNI_Cerebelum_9_R';       
%     ROI{24} = 'MNI_Cerebelum_Crus1_R';   
%     ROI{26} = 'MNI_Cerebelum_Crus2_R';   
%     ROI{28} = 'MNI_Cingulum_Ant_R';      
%     ROI{30} = 'MNI_Cingulum_Mid_R';      
%     ROI{32} = 'MNI_Cingulum_Post_R';           
%     ROI{58} = 'MNI_Hippocampus_R';  
%     ROI{72} = 'MNI_Pallidum_R'; 
%     ROI{108} = 'MNI_Thalamus_R';

    ROI{38} = 'MNI_Precentral_L'; 
    ROI{39} = 'MNI_Paracentral_Lobule_L';
    ROI{40} = 'MNI_Postcentral_L';
    ROI{41} = 'MNI_Frontal_Inf_Oper_L'; 
    ROI{42} = 'MNI_Frontal_Inf_Orb_L'; 
    ROI{43} = 'MNI_Frontal_Inf_Tri_L'; 
    ROI{44} = 'MNI_Frontal_Med_Orb_L';
    ROI{45} = 'MNI_Frontal_Mid_L';
    ROI{46} = 'MNI_Frontal_Sup_L'; 
    ROI{47} = 'MNI_Frontal_Mid_Orb_L'; 
    ROI{48} = 'MNI_Frontal_Sup_Medial_L';
    ROI{49} = 'MNI_Frontal_Sup_Orb_L'; 
    ROI{50} = 'MNI_Rolandic_Oper_L';  
    ROI{51} = 'MNI_Supp_Motor_Area_L'; 
    
    ROI{52} = 'MNI_Parietal_Inf_L'; 
    ROI{53} = 'MNI_Parietal_Sup_L'; 
    ROI{54} = 'MNI_Rectus_L'; 
    ROI{55} = 'MNI_ParaHippocampal_L'; 
    ROI{56} = 'MNI_SupraMarginal_L';
    
    ROI{57} = 'MNI_Precuneus_L'; 
    ROI{58} = 'MNI_Angular_L';   
    ROI{59} = 'MNI_Cuneus_L'; 
    ROI{60} = 'MNI_Calcarine_L'; 
    ROI{61} = 'MNI_Occipital_Inf_L'; 
    ROI{62} = 'MNI_Occipital_Mid_L'; 
    ROI{63} = 'MNI_Occipital_Sup_L';
    
    ROI{64} = 'MNI_Olfactory_L'; 
    ROI{65} = 'MNI_Lingual_L';
    ROI{66} = 'MNI_Heschl_L';
    ROI{67} = 'MNI_Fusiform_L'; 
    ROI{68} = 'MNI_Insula_L';
    ROI{69} = 'MNI_Temporal_Inf_L';
    ROI{70} = 'MNI_Temporal_Mid_L';
    ROI{71} = 'MNI_Temporal_Pole_Mid_L';
    ROI{72} = 'MNI_Temporal_Pole_Sup_L';
    ROI{73} = 'MNI_Temporal_Sup_L';
    ROI{74} = 'MNI_Putamen_L'; 
      
%     ROI{1} = 'MNI_Amygdala_L';      
%     ROI{7} = 'MNI_Caudate_L';  
%     ROI{9} = 'MNI_Cerebelum_10_L'; 
%     ROI{11} = 'MNI_Cerebelum_3_L'; 
%     ROI{13} = 'MNI_Cerebelum_4_5_L'; 
%     ROI{15} = 'MNI_Cerebelum_6_L'; 
%     ROI{17} = 'MNI_Cerebelum_7b_L'; 
%     ROI{19} = 'MNI_Cerebelum_8_L'; 
%     ROI{21} = 'MNI_Cerebelum_9_L';  
%     ROI{23} = 'MNI_Cerebelum_Crus1_L';   
%     ROI{25} = 'MNI_Cerebelum_Crus2_L'; 
%     ROI{27} = 'MNI_Cingulum_Ant_L';  
%     ROI{29} = 'MNI_Cingulum_Mid_L'; 
%     ROI{31} = 'MNI_Cingulum_Post_L';        
%     ROI{57} = 'MNI_Hippocampus_L'; 
%     ROI{71} = 'MNI_Pallidum_L';       
%     ROI{107} = 'MNI_Thalamus_L'; 

%     ROI{109} = 'MNI_Vermis_10';           
%     ROI{110} = 'MNI_Vermis_1_2';          
%     ROI{111} = 'MNI_Vermis_3';           
%     ROI{112} = 'MNI_Vermis_4_5';          
%     ROI{113} = 'MNI_Vermis_6';            
%     ROI{114} = 'MNI_Vermis_7';            
%     ROI{115} = 'MNI_Vermis_8';            
%     ROI{116} = 'MNI_Vermis_9'; 

end