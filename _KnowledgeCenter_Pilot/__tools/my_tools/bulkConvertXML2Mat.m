
all_xml_files = dir('*.xml');

nFiles = length(all_xml_files);


for iFile=1:nFiles
    
    try
    
        disp(strcat(int2str(iFile),':',all_xml_files(iFile).name));

        [xml_tree, ~, ~] = xml_read(all_xml_files(iFile).name);
    
        save(strcat(all_xml_files(iFile).name(1:end-4),'.mat'),'xml_tree');
        
    catch ME
        
        msg = getReport(ME);
        
        disp(strcat('ERROR:',int2str(iFile),':',all_xml_files(iFile).name));
        
        disp(msg);
    
    end
    
end

