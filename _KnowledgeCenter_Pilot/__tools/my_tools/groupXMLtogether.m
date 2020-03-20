
query = 'retina-1';

folder = 'D:\__PUBMED-references\utf-8\';

files = dir(strcat(folder,query,'\*.xml'));

for iFile=1:2
    
    [tree,~,~] = xml_read(strcat(folder,query,'\',files(iFile).name));
    
    all_tree(iFile).tree = tree;
    all_tree(iFile).fileName = files(iFile).name;
    
end