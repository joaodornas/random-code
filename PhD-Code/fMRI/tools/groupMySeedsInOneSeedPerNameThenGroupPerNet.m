function new_all_Nets_Densities = groupMySeedsInOneSeedPerNameThenGroupPerNet( Nets_Densities )

nNets = length(Nets_Densities.net);

for iNet=1:nNets
    
    new_all_Nets_Densities.net(iNet).Net_Seeds = cell(2,17);
    
end

for iNet=1:nNets
    
    nSeeds = size(Nets_Densities.net(iNet).Net_Seeds,1) - 1;
    
    new_all_Nets_Densities.net(iNet).Net_Seeds{1+1,1} = 'All_Seeds';
    
    new_data_numbers = Nets_Densities.net(iNet).Net_Seeds(1+1,2:2:16);
    
%     new_all_Nets_Densities.net(iNet).Net_Seeds(1+1,2:2:16);
    
    for iSeed=2:nSeeds
           
        new_data_numbers = Nets_Densities.net(iNet).Net_Seeds(iSeed+1,2:2:16);

        old_data_numbers = new_all_Nets_Densities.net(iNet).Net_Seeds(1+1,2:2:16);

        for i=1:length(new_data_numbers)

            if isempty(old_data_numbers{i}); old_number = 0; else old_number = old_data_numbers{i}; end
            if isempty(new_data_numbers{i}); new_number = 0; else new_number = new_data_numbers{i}; end

            sum_of_data = new_number + old_number;

            j = i*2;

            new_all_Nets_Densities.net(iNet).Net_Seeds{1+1,j} = sum_of_data;

        end

        for ii=[3 5 7 9 11 13 15 17]

            new_data = Nets_Densities.net(iNet).Net_Seeds{iSeed+1,ii};
            old_data = new_all_Nets_Densities.net(iNet).Net_Seeds{1+1,ii};

            new_all_Nets_Densities.net(iNet).Net_Seeds{1+1,ii} = [old_data;new_data];

        end
                
     end
            
end
    

for iNet=1:nNets
    
    for iaal=[3 5 7 9 11 13 15 17]

        aal_seeds = new_all_Nets_Densities.net(iNet).Net_Seeds{1+1,iaal};
        aal_summary = aal_per_seed_summary( aal_seeds );
        new_all_Nets_Densities.net(iNet).Net_Seeds{1+1,iaal} = aal_summary;

    end
    
end


end

