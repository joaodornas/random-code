

nK = size(moments_of_difficulty,1);
nMD = size(moments_of_difficulty,2);
nComb = size(moments_of_difficulty,3);
nColor = size(moments_of_difficulty,4);
nMOT = size(moments_of_difficulty,5);
nVelocities = size(moments_of_difficulty,6);

all_nMDColor = [];

velocities = [20 25 30 35 40 45 50];

for iVelocity=1:nVelocities

    this_MOT_nMDColor = [];
    
    for iMOT=1:nMOT
          
        results = zeros(nComb*nColor,6);

        iresults = 0;
        
        for iColor=1:nColor
            
            for iComb=1:nComb
            
                iresults = iresults + 1;
                
%                 nMDColor = sum(moments_of_difficulty(:,1,iComb,iColor,iMOT,iVelocity));
%                 nMDMotion = sum(moments_of_difficulty(:,2,iComb,iColor,iMOT,iVelocity));
%                 nMDArea = sum(moments_of_difficulty(:,3,iComb,iColor,iMOT,iVelocity));

                nMDColor = length(find(diff(moments_of_difficulty(:,1,iComb,iColor,iMOT,iVelocity))>0));
                nMDMotion = length(find(diff(moments_of_difficulty(:,2,iComb,iColor,iMOT,iVelocity))>0));
                nMDArea = length(find(diff(moments_of_difficulty(:,3,iComb,iColor,iMOT,iVelocity))>0));
            
                results(iresults,1) = iColor;
                results(iresults,2) = iComb;
                results(iresults,3) = nMDColor;
                results(iresults,4) = nMDMotion;
                results(iresults,5) = nMDArea;
                results(iresults,6) = nMDColor + nMDMotion + nMDArea;
                
            end
            
        end
        
        this_MOT_nMDColor = [this_MOT_nMDColor, nMDColor];
        
        MDMOTs.Velocity(iVelocity).MOT(iMOT).results = results;
        
    end
    
    [unique_nMDColor, nunique_nMDColor] = count_unique(this_MOT_nMDColor);
    
    all_nMDColor = [all_nMDColor; unique_nMDColor];
   
end

[unique_all_nMDColor, nunique_all_nMDColor] = count_unique(all_nMDColor);

idx_nunique_nMDColor = find(nunique_all_nMDColor == nVelocities);

max_common_nMDColor = max(unique_all_nMDColor(idx_nunique_nMDColor));

all_nMDSum = [];

for iVelocity=1:nVelocities

    this_MOT_nMDSum = [];
    
    for iMOT=1:nMOT

        iresults = 0;
        
        for iColor=1:nColor
            
            for iComb=1:nComb
            
                iresults = iresults + 1;
                
                if MDMOTs.Velocity(iVelocity).MOT(iMOT).results(iresults,3) == max_common_nMDColor
                    
                    this_MOT_nMDSum = [this_MOT_nMDSum, MDMOTs.Velocity(iVelocity).MOT(iMOT).results(iresults,6)];
                    
                end
                
            end
            
        end
        
    end
    
    [unique_nMDSum, nunique_nMDSum] = count_unique(this_MOT_nMDSum);
    
    all_nMDSum = [all_nMDSum; unique_nMDSum];
    
end

[unique_all_nMDSum, nunique_all_nMDSum] = count_unique(all_nMDSum);

idx_nunique_nMDSum = find(nunique_all_nMDSum == nVelocities);

max_common_nMDSum = max(unique_all_nMDSum(idx_nunique_nMDSum));

TotalMDs = zeros(nComb*2,nVelocities);

for iVelocity=1:nVelocities

    for iMOT=1:nMOT
              
        if max_common_nMDColor == MDMOTs.Velocity(iVelocity).MOT(iMOT).results(1,3)
            
            total = MDMOTs.Velocity(iVelocity).MOT(iMOT).results(:,6);
            idx_max_sum = find(total == max_common_nMDSum);
            
            if ~isempty(idx_max_sum)
                
                color = MDMOTs.Velocity(iVelocity).MOT(iMOT).results(idx_max_sum(1),1);
                comb = MDMOTs.Velocity(iVelocity).MOT(iMOT).results(idx_max_sum(1),2);
            
                MDTargetsCombinations{iVelocity,1} = color;
                MDTargetsCombinations{iVelocity,2} = comb;
                
                MDTargetsCombinations{iVelocity,3} = iMOT;
                
                prefix = 'MOT-Analyze-V2-11-m-';
                analyzeareadata = load(strcat(prefix,'Area','-',int2str(velocities(iVelocity)),'-',int2str(iMOT),'.mat'));
              
                allColorComb = analyzeareadata.allColorComb;
                
                MDTargetsCombinations{iVelocity,4} = allColorComb(color).comb(comb,:);
            
            end
            
            TotalMDs(:,iVelocity) = MDMOTs.Velocity(iVelocity).MOT(iMOT).results(:,6);
            
        end
        
    end
    
end

for iComb=1:nComb*2
    
    MSQDIFFMOD(iComb) = sqrt( (TotalMDs(iComb,1) - TotalMDs(iComb,2))^2 + (TotalMDs(iComb,1) - TotalMDs(iComb,3))^2 + (TotalMDs(iComb,1) - TotalMDs(iComb,4))^2 + (TotalMDs(iComb,1) - TotalMDs(iComb,5))^2 + (TotalMDs(iComb,1) - TotalMDs(iComb,6))^2 + (TotalMDs(iComb,1) - TotalMDs(iComb,7))^2 );
                 
end

idx_comb_msqdiff = find(MSQDIFFMOD == min(MSQDIFFMOD));



