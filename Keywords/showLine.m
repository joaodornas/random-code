

xmlfile = input('XML File:','s');

whichLine = input('Line:');

xml = fopen(xmlfile, 'r', 'l', 'UTF-8');

l = 0;

while l <= whichLine
    
    line = fgetl(xml);
    l = l + 1;
    
    if l == whichLine
        
        disp(line);
        
    end
    
end

fclose('all');