

xmlfile = input('XML File:','s');

xml = fopen(xmlfile, 'r', 'l', 'UTF-8');

fileToWrite = fopen(strcat(xmlfile(1:end-4),'-corrected.xml'), 'w', 'l', 'UTF-8');

xml_string = '<?xml version="1.0"?>';
doc_string = '<!DOCTYPE ';
start_of_section_string = '<PubmedArticleSet>';
end_of_section_string = '</PubmedArticleSet>';

l = 0;

while ~feof(xml)
    
    line = fgetl(xml);
    l = l + 1;
    
    if l <= 3
        
        fwrite(fileToWrite, line, 'char', 0, 's');
        fprintf(fileToWrite, '\n');

    else
    
        idx_xml_string = [];
        idx_doc_string = [];
        idx_start_of_section_string = [];
        idx_end_of_section_string = [];

        idx_test_strings = [];

        idx_xml_string = strfind(line,xml_string);
        idx_doc_string = strfind(line,doc_string);
        idx_start_of_section_string = strfind(line,start_of_section_string);
        idx_end_of_section_string = strfind(line,end_of_section_string);

        testStrings = [~isempty(idx_xml_string) ~isempty(idx_doc_string) ~isempty(idx_start_of_section_string) ~isempty(idx_end_of_section_string)];

        idx_test_strings = find(testStrings);

        if isempty(idx_test_strings)
            
            line = strrep(line,'?',' '); 

            fwrite(fileToWrite, line, 'char', 0, 's');

            fprintf(fileToWrite, '\n');

        end
    
    end
 
end

fwrite(fileToWrite, end_of_section_string, 'char', 0, 's');
fprintf(fileToWrite, '\n');

fclose(fileToWrite);

fclose('all');