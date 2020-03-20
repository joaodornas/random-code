function lowhigh_FC_ICA_per_ROI


settings_jan_0805;
%settings_elena_2905;

%doTheMath(settings);

plotResults(settings);

%plotConnections(settings);

end

function doTheMath(settings)

preprocessed_folder = strcat(settings.folders.main,'\',settings.folders.experiment,'\',settings.folders.subject,'\',settings.folders.preprocessed);

ICA_folder{1} = strcat(preprocessed_folder,'\',settings.functional.mot4.folder.main,'\',settings.functional.mot4.run(1).folder,'\','FSL','\','ICA');
ICA_folder{2} = strcat(preprocessed_folder,'\',settings.functional.mot4.folder.main,'\',settings.functional.mot4.run(2).folder,'\','FSL','\','ICA');

ICA_folder{3} = strcat(preprocessed_folder,'\',settings.functional.mot2.folder.main,'\',settings.functional.mot2.run(1).folder,'\','FSL','\','ICA');
ICA_folder{4} = strcat(preprocessed_folder,'\',settings.functional.mot2.folder.main,'\',settings.functional.mot2.run(2).folder,'\','FSL','\','ICA');

ICA_folder{5} = strcat(preprocessed_folder,'\',settings.functional.restingstate.folder.main,'\',settings.functional.restingstate.run(1).folder,'\','FSL','\','ICA');
ICA_folder{6} = strcat(preprocessed_folder,'\',settings.functional.restingstate.folder.main,'\',settings.functional.restingstate.run(2).folder,'\','FSL','\','ICA');

ICA_label{1} = 'High-Attention-Run-1';
ICA_label{2} = 'High-Attention-Run-2';
ICA_label{3} = 'Low-Attention-Run-1';
ICA_label{4} = 'Low-Attention-Run-2';
ICA_label{5} = 'RestingState-Run-1';
ICA_label{6} = 'RestingState-Run-2';

nTR = 331;
nFalse_ic = 10;

for ifolder=1:length(ICA_folder)
    
    disp(ICA_label{ifolder});
    
    nNodes = 90;
    
    allComponents = [];
    allPercentComponents = [];
    for iNode=1:nNodes
        
        disp(strcat('ROI:',int2str(iNode)));
        
        if iNode ~= 87
            
            ROI(iNode).ic = load(strcat(ICA_folder{ifolder},'\','aal2std_',int2str(iNode),'\','melodic_mix'));
        
        else
            
            ROI(iNode).ic = zeros(nTR,nFalse_ic);
            
        end
        
        ROI(iNode).nIC = size(ROI(iNode).ic,2);
        
        allComponents = [allComponents,ROI(iNode).ic];
        
        if iNode ~= 87
            
            ic_var = textread(strcat(ICA_folder{ifolder},'\','aal2std_',int2str(iNode),'\','eigenvalues_percent'),'%s');
            ic_var = str2double(ic_var);

            variance_explained = 0.5;
            ICA_n(ifolder).variance_explained = variance_explained;
            explained = 0;
            nPercentComponents = 0;
            while explained < variance_explained

                nPercentComponents = nPercentComponents + 1;

                explained = ic_var(nPercentComponents);

            end

        else
            
            nPercentComponents = nFalse_ic;
            
        end
        
        ICA_n(ifolder).nIC(iNode) = size(ROI(iNode).ic,2);
        ICA_n(ifolder).nICpercent(iNode) = nPercentComponents;

        allPercentComponents = [allPercentComponents,ROI(iNode).ic(:,1:nPercentComponents)];
        
    end
    
    disp('doing rho');
    
    [ICA_FC(ifolder).rho,ICA_FC(ifolder).pval] = corr(allComponents);
    [ICA_FC(ifolder).percent_rho,ICA_FC(ifolder).percent_pval] = corr(allPercentComponents);
 
end

save(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-per-ROI-All-Runs','.mat'),'ICA_label','ICA_folder','ICA_n','ICA_FC','-v7.3');

end

function plotResults(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-per-ROI-All-Runs','.mat'));

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

nNodes = 90;

for ifolder=1:length(ICA_folder)
    
   f = figure;
   
   imagesc(ICA_FC(ifolder).rho,[min_all max_all]);
   
   xlabel('AAL Regions ICs');
   ylabel('AAL Regions ICs');
   
   title(ICA_label{ifolder});
   
   colorbar;
   
%    hold on;
%     
%    for iNode=1:nNodes
%         
%         p1(1) = sum(ICA_n(ifolder).nIC(1:iNode));
%         p2(1) = 1;
%         p1(2) = sum(ICA_n(ifolder).nIC(1:iNode));
%         p2(2) = size(ICA_FC(ifolder).rho,1);
%    
%         plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','k','LineWidth',2);
%     
%         p1(1) = 1;
%         p2(1) = sum(ICA_n(ifolder).nIC(1:iNode));
%         p1(2) = size(ICA_FC(ifolder).rho,1);
%         p2(2) = sum(ICA_n(ifolder).nIC(1:iNode));
%    
%         plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','k','LineWidth',2);
%    
%    end
   
   print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-per-ROI','-',ICA_label{ifolder},'.jpeg'));
   print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-per-ROI','-',ICA_label{ifolder},'.eps'));
   print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-per-ROI','-',ICA_label{ifolder},'.pdf'));
   
   
   f = figure;
   
   imagesc(ICA_FC(ifolder).percent_rho,[min_percent_all max_percent_all]);
   
   xlabel('AAL Regions ICs');
   ylabel('AAL Regions ICs');
   
   title(strcat(ICA_label{ifolder},':',int2str(round(ICA_n(ifolder).variance_explained*100)),'% of variance'));
   
   colorbar;
   
%    hold on;
%    
%    for iNode=1:nNodes
%        
%         p1(1) = sum(ICA_n(ifolder).nICpercent(1:iNode));
%         p2(1) = 1;
%         p1(2) = sum(ICA_n(ifolder).nICpercent(1:iNode));
%         p2(2) = size(ICA_FC(ifolder).percent_rho,1);
%    
%         plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','k','LineWidth',2);
%     
%         p1(1) = 1;
%         p2(1) = sum(ICA_n(ifolder).nICpercent(1:iNode));
%         p1(2) = size(ICA_FC(ifolder).percent_rho,1);
%         p2(2) = sum(ICA_n(ifolder).nICpercent(1:iNode));
%    
%         plot([p1(1),p1(2)],[p2(1),p2(2)],'Color','k','LineWidth',2);
%    
%    end
   
   print(f,'-djpeg',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-per-ROI-percent','-',ICA_label{ifolder},'.jpeg'));
   print(f,'-depsc',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-per-ROI-percent','-',ICA_label{ifolder},'.eps'));
   print(f,'-dpdf',strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-per-ROI-percent','-',ICA_label{ifolder},'.pdf'));
   
   close all,
   
end

end

function plotConnections(settings)

load(strcat(settings.folders.experiment,'-',settings.folders.subject,'-','ICA-per-ROI-All-Runs','.mat'));

nFolders = 6;
nNodes = 90;

for ifolder=1:nFolders
   
   for iNode=1:nNodes
       
      ICA_n(ifolder).nICpercent(iNode);
       
       
   end
    
    
end

end

