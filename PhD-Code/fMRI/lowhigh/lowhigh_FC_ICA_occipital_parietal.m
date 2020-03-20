function lowhigh_FC_ICA_occipital_parietal


settings_jan_0805;
%settings_elena_2905;

%doTheMath(settings);

plotResults(settings);

end

function doTheMath(settings)

preprocessed_folder = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',settings.folders.preprocessed);

ICA_folder{1} = strcat(preprocessed_folder,'\',settings.functional.mot4.folder.main,'\',settings.functional.mot4.run(1).folder,'\','FSL','\','ICA_area');
ICA_folder{2} = strcat(preprocessed_folder,'\',settings.functional.mot4.folder.main,'\',settings.functional.mot4.run(2).folder,'\','FSL','\','ICA_area');

ICA_folder{3} = strcat(preprocessed_folder,'\',settings.functional.mot2.folder.main,'\',settings.functional.mot2.run(1).folder,'\','FSL','\','ICA_area');
ICA_folder{4} = strcat(preprocessed_folder,'\',settings.functional.mot2.folder.main,'\',settings.functional.mot2.run(2).folder,'\','FSL','\','ICA_area');

ICA_folder{5} = strcat(preprocessed_folder,'\',settings.functional.restingstate.folder.main,'\',settings.functional.restingstate.run(1).folder,'\','FSL','\','ICA_area');
ICA_folder{6} = strcat(preprocessed_folder,'\',settings.functional.restingstate.folder.main,'\',settings.functional.restingstate.run(2).folder,'\','FSL','\','ICA_area');

ICA_label{1} = 'High-Attention-Run-1';
ICA_label{2} = 'High-Attention-Run-2';
ICA_label{3} = 'Low-Attention-Run-1';
ICA_label{4} = 'Low-Attention-Run-2';
ICA_label{5} = 'RestingState-Run-1';
ICA_label{6} = 'RestingState-Run-2';

for ifolder=1:length(ICA_folder)
    
    disp(ICA_label{ifolder});
    
    occipital_ic = load(strcat(ICA_folder{ifolder},'\','occipital','\','melodic_mix'));
    parietal_ic = load(strcat(ICA_folder{ifolder},'\','parietal','\','melodic_mix'));
    
    ICA_n(ifolder).nOccipital = size(occipital_ic,2);
    ICA_n(ifolder).nParietal = size(parietal_ic,2);
    
    allComponents = [occipital_ic,parietal_ic];
    
    [ICA_FC(ifolder).rho,ICA_FC(ifolder).pval] = corr(allComponents);
 
    occipital_var = textread(strcat(ICA_folder{ifolder},'\','occipital','\','eigenvalues_percent'),'%s');
    parietal_var = textread(strcat(ICA_folder{ifolder},'\','parietal','\','eigenvalues_percent'),'%s');
    
    occipital_var = str2double(occipital_var);
    parietal_var = str2double(parietal_var);
    
    variance_explained = 0.5;
    ICA_n(ifolder).variance_explained = variance_explained;
    explained = 0;
    ic = 0;
    while explained < variance_explained
        
        ic = ic + 1;
        
        explained = occipital_var(ic);
        
    end
    
    ICA_n(ifolder).nICpercentOccipital = ic;
    
    variance_explained = 0.5;
    explained = 0;
    ic = 0;
    while explained < variance_explained
        
        ic = ic + 1;
        
        explained = parietal_var(ic);
        
    end
    
    ICA_n(ifolder).nICpercentParietal = ic;
    
    percentOccipital_ic = occipital_ic(:,1:ICA_n(ifolder).nICpercentOccipital);
    percentParietal_ic = parietal_ic(:,1:ICA_n(ifolder).nICpercentParietal);
    
    percentComponents = [percentOccipital_ic,percentParietal_ic];
    
    [ICA_FC(ifolder).percent_rho,ICA_FC(ifolder).percent_pval] = corr(percentComponents);
 
end

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal-All-Runs','.mat'),'ICA_label','ICA_folder','ICA_n','ICA_FC');

end


function plotResults(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal-All-Runs','.mat'));

for ifolder=1:length(ICA_folder)
    
   max_rho(ifolder) = max(max(ICA_FC(ifolder).rho));
   min_rho(ifolder) = min(min(ICA_FC(ifolder).rho));
   
   max_percent_rho(ifolder) =  max(max(ICA_FC(ifolder).percent_rho));
   min_percent_rho(ifolder) = min(min(ICA_FC(ifolder).percent_rho));
    
end

max_all = max(max_rho);
min_all = min(min_rho);

max_percent_all = min(max_percent_rho);
min_percent_all = min(min_percent_rho);

for ifolder=1:length(ICA_folder)
    
   f = figure;
   
   imagesc(ICA_FC(ifolder).rho,[min_all max_all]);
   
   xlabel('occipital - parietal');
   ylabel('parietal - occipital');
   
   title(ICA_label{ifolder});
   
   colorbar;
   
   hold on;
   
   p1(1) = ICA_n(ifolder).nOccipital;
   p2(1) = 1;
   p1(2) = ICA_n(ifolder).nOccipital;
   p2(2) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   p1(1) = 1;
   p2(1) = ICA_n(ifolder).nOccipital;
   p1(2) = ICA_n(ifolder).nOccipital + ICA_n(ifolder).nParietal;
   p2(2) = ICA_n(ifolder).nOccipital;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal','-',ICA_label{ifolder},'.jpeg'));
   print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal','-',ICA_label{ifolder},'.eps'));
   print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal','-',ICA_label{ifolder},'.pdf'));
   
   
   f = figure;
   
   imagesc(ICA_FC(ifolder).percent_rho,[min_percent_all max_percent_all]);
   
   xlabel('occipital - parietal');
   ylabel('parietal - occipital');
   
   title(strcat(ICA_label{ifolder},':',int2str(round(ICA_n(ifolder).variance_explained*100)),'% of variance'));
   
   colorbar;
   
   hold on;
   
   p1(1) = ICA_n(ifolder).nICpercentOccipital;
   p2(1) = 1;
   p1(2) = ICA_n(ifolder).nICpercentOccipital;
   p2(2) = ICA_n(ifolder).nICpercentOccipital + ICA_n(ifolder).nICpercentParietal;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   p1(1) = 1;
   p2(1) = ICA_n(ifolder).nICpercentOccipital;
   p1(2) = ICA_n(ifolder).nICpercentOccipital + ICA_n(ifolder).nICpercentParietal;
   p2(2) = ICA_n(ifolder).nICpercentOccipital;
   
   plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','r','LineWidth',2);
   
   print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal-percent','-',ICA_label{ifolder},'.jpeg'));
   print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal-percent','-',ICA_label{ifolder},'.eps'));
   print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-Occipital-Parietal-percent','-',ICA_label{ifolder},'.pdf'));
   
end

end

