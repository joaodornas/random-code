function aal_summary = aal_per_seed_summary( aal_seeds )

    aal_names = aal_seeds(:,1);
    
    idx = [];
    for i=1:length(aal_names)
        if isempty(aal_names{i})
        idx = [idx,i];
        end
    end
       
    aal_names(idx) = [];
    aal_seeds(idx,:) = [];
    
    [C,IA,IC] = unique(aal_names);
    
    voxels = zeros(1,length(C));
    
    for iC=1:length(C)
        
        for iAAL=1:size(aal_seeds)
            
            if strcmp(C(iC),aal_seeds(iAAL,1))
                
                voxels(iC) = voxels(iC) + aal_seeds{iAAL,2};
                
            end
            
        end
        
    end
    
    aal_summary = cell(length(C),2);
    
    for iC=1:length(C)
        
        aal_summary{iC,1} = C{iC};
        aal_summary{iC,2} = voxels(iC);
        
    end

end

