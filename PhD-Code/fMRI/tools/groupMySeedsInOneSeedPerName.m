function new_Nets_Densities = groupMySeedsInOneSeedPerName( Nets_Densities )

nNets = length(Nets_Densities.net);

for iNet=1:nNets
    
    new_Nets_Densities.net(iNet).StatsResults = Nets_Densities.net(iNet).StatsResults;
    new_Nets_Densities.net(iNet).network_label = Nets_Densities.net(iNet).network_label;
   
    nSeeds = size(Nets_Densities.net(iNet).Net_Seeds,1) - 1;
    
    all_seeds_names = Nets_Densities.net(iNet).Net_Seeds(2:end,1);
  
    for iSeed=1:nSeeds
       
        name = all_seeds_names{iSeed};
     
        new_name = cleanSeedName(name);
        
        all_seeds_names{iSeed} = new_name;
        
    end
    
    [net(iNet).unique_seeds, net(iNet).count_seeds] = count_unique(all_seeds_names);

end

for iNet=1:nNets
    
    net(iNet).nGroups = length(net(iNet).count_seeds);
    net(iNet).first_time_group = zeros(1,net(iNet).nGroups);
    
    new_Nets_Densities.net(iNet).Net_Seeds = cell(net(iNet).nGroups+1,size(Nets_Densities.net(iNet).Net_Seeds,2));
    new_Nets_Densities.net(iNet).Net_Seeds(1,:) = Nets_Densities.net(iNet).Net_Seeds(1,:);
    
end

for iNet=1:nNets
    
    nSeeds = size(Nets_Densities.net(iNet).Net_Seeds,1) - 1;
    
    for iSeed=1:nSeeds
        
        seed_info = Nets_Densities.net(iNet).Net_Seeds(iSeed+1,:);
        
        seed_name = seed_info(1);
        
        cleaned_seed_name = cleanSeedName(seed_name{1});
        
        for iGroup=1:net(iNet).nGroups
            
            if ~isempty(strfind(cleaned_seed_name,net(iNet).unique_seeds{iGroup}))
                
                iiSeed = iGroup;
                
                new_Nets_Densities.net(iNet).Net_Seeds{iiSeed+1,1} = cleaned_seed_name;
                   
                new_data_numbers = Nets_Densities.net(iNet).Net_Seeds(iSeed+1,2:2:16);

                old_data_numbers = new_Nets_Densities.net(iNet).Net_Seeds(iiSeed+1,2:2:16);

                for i=1:length(new_data_numbers)
                    
                    if isempty(old_data_numbers{i}); old_number = 0; else old_number = old_data_numbers{i}; end
                    if isempty(new_data_numbers{i}); new_number = 0; else new_number = new_data_numbers{i}; end

                    sum_of_data = new_number + old_number;

                    j = i*2;

                    new_Nets_Densities.net(iNet).Net_Seeds{iiSeed+1,j} = sum_of_data;

                end

                for ii=[3 5 7 9 11 13 15 17]

                    new_data = Nets_Densities.net(iNet).Net_Seeds{iSeed+1,ii};
                    old_data = new_Nets_Densities.net(iNet).Net_Seeds{iiSeed+1,ii};

                    new_Nets_Densities.net(iNet).Net_Seeds{iiSeed+1,ii} = [old_data;new_data];

                end
                
            end
            
        end
    
    end
    
end

for iNet=1:nNets
    
    nSeeds = size(new_Nets_Densities.net(iNet).Net_Seeds,1) - 1;
    
    for iSeed=1:nSeeds
        
        for iaal=[3 5 7 9 11 13 15 17]

            aal_seeds = new_Nets_Densities.net(iNet).Net_Seeds{iSeed+1,iaal};
            aal_summary = aal_per_seed_summary( aal_seeds );
            new_Nets_Densities.net(iNet).Net_Seeds{iSeed+1,iaal} = aal_summary;

        end
    
    end
    
end
 
  



end

