
all_xml_files = dir('*.xml');

nFiles = length(all_xml_files);


for iFile=1:nFiles
    
    try
    
        disp(strcat(int2str(iFile),':',all_xml_files(iFile).name));

        fin = fopen(all_xml_files(iFile).name,'r');
        fout = fopen(strcat(all_xml_files(iFile).name,'-tmp'),'w');
        
        while ~feof(fin)
            
            s = fgetl(fin);
            s = strrep(s,' & ',' &amp; ');
            s = strrep(s,' < ',' &lt; ');
            s = strrep(s,' > ',' &gt; ');
            s = s(s>31);
            
            fprintf(fout,'%s',s);
            
        end
        
        fclose(fin);
        fclose(fout);
        
        system(sprintf('move "%s" "%s"',strcat(all_xml_files(iFile).name,'-tmp'),all_xml_files(iFile).name));
        
    catch ME
        
        msg = getReport(ME);
        
        disp(strcat('ERROR:',int2str(iFile),':',all_xml_files(iFile).name));
        
        disp(msg);
    
    end
    
end