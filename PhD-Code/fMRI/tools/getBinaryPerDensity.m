function binary = getBinaryPerDensity(Nets_Densities,binary,density_ID,iSeed,iNet,nSeeds,nAAL_ROI,AAL_img,AAL_ROI)

    last_seeds = 0;
    for iiNet=1:(iNet-1)
        
        last_seeds = size(Nets_Densities.net(iiNet).Net_Seeds,1) - 1 + last_seeds;
        
    end
        
    if Nets_Densities.net(iNet).Net_Seeds{iSeed+1,density_ID} > 0;
            
            aal_cell = Nets_Densities.net(iNet).Net_Seeds{iSeed+1,density_ID+1};
            
            for i=2:size(aal_cell,1)
                
                name = aal_cell{i,1};
                
                for iROI=1:nAAL_ROI
                    
                   aal_name = AAL_ROI(iROI).Nom_L;
                   
                   if strcmp(name,aal_name)
                       
                       binary(iROI,nAAL_ROI+last_seeds+iSeed) = 1;
                       binary(nAAL_ROI+last_seeds+iSeed,iROI) = 1;
                       
                   end
                   
                end
                
            end
            
     end
        
end