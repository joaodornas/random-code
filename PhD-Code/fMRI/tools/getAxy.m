function [neg_up, neg_down, pos_up, pos_down, xy] = getAxy(Nets_Densities,pnu,pnd,ppu,ppd,anu,and,apu,apd,x,y,nSeeds,nAAL_ROI,Net)

for i=1:length(x)
    
      xy(i,:) = [x(i) y(i)];
%     xy(i,1) = x(i);
%     xy(i,2) = y(i);
    
end

pcriterion = 0.01;

%%% REMOVE NON-SIGNIFICANT BASED ON WILCOXON

all_networks_labels = {'DAN','VAN','SMN','VIS','FPC','LAN','AUD','DMN'};

nNet = 8;

neg_up = anu;
neg_down = and;
pos_up = apu;
pos_down = apd;

for iNet=1:nNet
    
    nNet_Seeds = Net(iNet).nSeeds;
    
    last_seeds = 0;
    for iiNet=1:(iNet-1)
        
        last_seeds = size(Nets_Densities.net(iiNet).Net_Seeds,1) - 1 + last_seeds;
        
    end
    
    if Nets_Densities.net(iNet).StatsResults{2,3} > pcriterion
        
        neg_down((nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds),:) = 0;
        neg_down(:,(nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds)) = 0;
        
    end
        
    if Nets_Densities.net(iNet).StatsResults{3,3} > pcriterion
        
        neg_up((nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds),:) = 0;
        neg_up(:,(nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds)) = 0;
        
    end
        
    if Nets_Densities.net(iNet).StatsResults{4,3} > pcriterion
        
        pos_down((nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds),:) = 0;
        pos_down(:,(nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds)) = 0;
        
    end
    
    if Nets_Densities.net(iNet).StatsResults{5,3} > pcriterion
        
        pos_up((nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds),:) = 0;
        pos_up(:,(nAAL_ROI+last_seeds+1):(nAAL_ROI+last_seeds+nNet_Seeds)) = 0;
        
    end
  
end

end


