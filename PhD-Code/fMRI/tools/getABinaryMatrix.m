function [pnu,pnd,ppu,ppd,anu,and,apu,apd] =  getABinaryMatrix(Nets_Densities,nSeeds,nAAL_ROI,AAL_img,AAL_ROI)

nNets = 8;

pnu = zeros(nAAL_ROI+nSeeds);
pnd = zeros(nAAL_ROI+nSeeds);
ppu = zeros(nAAL_ROI+nSeeds);
ppd = zeros(nAAL_ROI+nSeeds);

anu = zeros(nAAL_ROI+nSeeds);
and = zeros(nAAL_ROI+nSeeds);
apu = zeros(nAAL_ROI+nSeeds);
apd = zeros(nAAL_ROI+nSeeds);

for iNet=1:nNets
    
    n_net_Seeds = size(Nets_Densities.net(iNet).Net_Seeds,1) - 1;

    for iSeed=1:n_net_Seeds
          
        pnu = getBinaryPerDensity(Nets_Densities,pnu,4,iSeed,iNet,nSeeds,nAAL_ROI,AAL_img,AAL_ROI);
        pnd = getBinaryPerDensity(Nets_Densities,pnd,2,iSeed,iNet,nSeeds,nAAL_ROI,AAL_img,AAL_ROI);
        ppu = getBinaryPerDensity(Nets_Densities,ppu,8,iSeed,iNet,nSeeds,nAAL_ROI,AAL_img,AAL_ROI);
        ppd = getBinaryPerDensity(Nets_Densities,ppd,6,iSeed,iNet,nSeeds,nAAL_ROI,AAL_img,AAL_ROI);
        
        anu = getBinaryPerDensity(Nets_Densities,anu,12,iSeed,iNet,nSeeds,nAAL_ROI,AAL_img,AAL_ROI);
        and = getBinaryPerDensity(Nets_Densities,and,10,iSeed,iNet,nSeeds,nAAL_ROI,AAL_img,AAL_ROI);
        apu = getBinaryPerDensity(Nets_Densities,apu,16,iSeed,iNet,nSeeds,nAAL_ROI,AAL_img,AAL_ROI);
        apd = getBinaryPerDensity(Nets_Densities,apd,14,iSeed,iNet,nSeeds,nAAL_ROI,AAL_img,AAL_ROI);
        
    end
  
end

end