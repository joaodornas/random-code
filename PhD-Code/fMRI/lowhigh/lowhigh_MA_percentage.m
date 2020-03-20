
nos = load('AAL-nodes.mat');
nodes = nos.nodes;

%percentage = 20/100;
percentage = 45/100;

nNodes = size(nodes,1);

MOT4RestingState_MA_edges = zeros(nNodes,nNodes);
MOT2RestingState_MA_edges = zeros(nNodes,nNodes);
MOT4MOT2_MA_edges = zeros(nNodes,nNodes);

for iStructure=1:nNodes
    
    for iiStructure=iStructure:nNodes
        
        MOT4_MA_value = MOT4_MA{iStructure+1,iiStructure+1};
        MOT2_MA_value = MOT2_MA{iStructure+1,iiStructure+1};
        RestingState_MA_value = RestingState_MA{iStructure+1,iiStructure+1};
        
        if abs(MOT4_MA_value - RestingState_MA_value) > percentage
            
            MOT4RestingState_MA_edges(iStructure,iiStructure) = 1*sign(MOT4_MA_value - RestingState_MA_value);
            MOT4RestingState_MA_edges(iiStructure,iStructure) = 1*sign(MOT4_MA_value - RestingState_MA_value);
        
        end
        
        if abs(MOT2_MA_value - RestingState_MA_value) > percentage
            
            MOT2RestingState_MA_edges(iStructure,iiStructure) = 1*sign(MOT2_MA_value - RestingState_MA_value);
            MOT2RestingState_MA_edges(iiStructure,iStructure) = 1*sign(MOT2_MA_value - RestingState_MA_value);
        
        end
           
        if abs(MOT4_MA_value - MOT2_MA_value) > percentage
            
            MOT4MOT2_MA_edges(iStructure,iiStructure) = 1*sign(MOT4_MA_value - MOT2_MA_value);
            MOT4MOT2_MA_edges(iiStructure,iStructure) = 1*sign(MOT4_MA_value - MOT2_MA_value);
        
        end
           
    end
    
end

MOT4RestingState_nodes = cell.empty;
MOT2RestingState_nodes = cell.empty;
MOT4MOT2_nodes = cell.empty;

idx_MOT4RestingState = [];
idx_MOT2RestingState = [];
idx_MOT4MOT2 = [];

for iStructure=1:nNodes
    
    line = MOT4RestingState_MA_edges(iStructure,:);
    
    if sum(line) == 0
        
        idx_MOT4RestingState = [idx_MOT4RestingState, iStructure];
        
    end
    
    line = MOT2RestingState_MA_edges(iStructure,:);
    
    if sum(line) == 0
        
        idx_MOT2RestingState = [idx_MOT2RestingState, iStructure];
        
    end
    
    line = MOT4MOT2_MA_edges(iStructure,:);
    
    if sum(line) == 0
        
        idx_MOT4MOT2 = [idx_MOT4MOT2, iStructure];
        
    end
    
end

MOT4RestingState_nodes = nodes;
MOT2RestingState_nodes = nodes;
MOT4MOT2_nodes = nodes;

if ~isempty(idx_MOT4RestingState)
    
    MOT4RestingState_nodes(idx_MOT4RestingState,:) = [];
    
    MOT4RestingState_MA_edges(idx_MOT4RestingState,:) = [];
    MOT4RestingState_MA_edges(:,idx_MOT4RestingState) = [];
    
end

if ~isempty(idx_MOT2RestingState)
    
    MOT2RestingState_nodes(idx_MOT2RestingState,:) = [];
    
    MOT2RestingState_MA_edges(idx_MOT2RestingState,:) = [];
    MOT2RestingState_MA_edges(:,idx_MOT2RestingState) = [];
    
end

if ~isempty(idx_MOT4MOT2)
    
    MOT4MOT2_nodes(idx_MOT4MOT2,:) = [];
    
    MOT4MOT2_MA_edges(idx_MOT4MOT2,:) = [];
    MOT4MOT2_MA_edges(:,idx_MOT4MOT2) = [];
    
end

MOT4RestingState_nodes_file = fopen(strcat('Low-High-',settings.folders.subject,'-MOT4RestingState-MA-voxels',num2str(percentage),'-nodes.txt'),'wt','l','UTF-8');

for iStructure=1:size(MOT4RestingState_nodes,1)
    
    fprintf(MOT4RestingState_nodes_file,strcat(MOT4RestingState_nodes{iStructure,1},'%s',MOT4RestingState_nodes{iStructure,2},'%s',MOT4RestingState_nodes{iStructure,3},'%s',MOT4RestingState_nodes{iStructure,4},'%s',MOT4RestingState_nodes{iStructure,5},'%s',MOT4RestingState_nodes{iStructure,6},'\n'),' ',' ',' ',' ',' ');
    
end

fclose(MOT4RestingState_nodes_file);

MOT2RestingState_nodes_file = fopen(strcat('Low-High-',settings.folders.subject,'-MOT2RestingState-MA-voxels',num2str(percentage),'-nodes.txt'),'wt','l','UTF-8');

for iStructure=1:size(MOT2RestingState_nodes,1)
    
    fprintf(MOT2RestingState_nodes_file,strcat(MOT2RestingState_nodes{iStructure,1},'%s',MOT2RestingState_nodes{iStructure,2},'%s',MOT2RestingState_nodes{iStructure,3},'%s',MOT2RestingState_nodes{iStructure,4},'%s',MOT2RestingState_nodes{iStructure,5},'%s',MOT2RestingState_nodes{iStructure,6},'\n'),' ',' ',' ',' ',' ');
    
end

fclose(MOT2RestingState_nodes_file);

MOT4MOT2_nodes_file = fopen(strcat('Low-High-',settings.folders.subject,'-MOT4MOT2-MA-voxels',num2str(percentage),'-nodes.txt'),'wt','l','UTF-8');

for iStructure=1:size(MOT4MOT2_nodes,1)
    
    fprintf(MOT4MOT2_nodes_file,strcat(MOT4MOT2_nodes{iStructure,1},'%s',MOT4MOT2_nodes{iStructure,2},'%s',MOT4MOT2_nodes{iStructure,3},'%s',MOT4MOT2_nodes{iStructure,4},'%s',MOT4MOT2_nodes{iStructure,5},'%s',MOT4MOT2_nodes{iStructure,6},'\n'),' ',' ',' ',' ',' ');
    
end

fclose(MOT4MOT2_nodes_file);

dlmwrite(strcat('Low-High-',settings.folders.subject,'-MOT4RestingState-MA-voxels',num2str(percentage),'-edges.txt'),MOT4RestingState_MA_edges,'-append','delimiter',' ','newline','pc');
dlmwrite(strcat('Low-High-',settings.folders.subject,'-MOT2RestingState-MA-voxels',num2str(percentage),'-edges.txt'),MOT2RestingState_MA_edges,'-append','delimiter',' ','newline','pc');
dlmwrite(strcat('Low-High-',settings.folders.subject,'-MOT4MOT2-MA-voxels',num2str(percentage),'-edges.txt'),MOT4MOT2_MA_edges,'-append','delimiter',' ','newline','pc');
