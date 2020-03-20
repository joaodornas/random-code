function countLines(xmlfile)

xml = fopen(xmlfile, 'r', 'l', 'UTF-8');

l = 0;
threshold = 100000;
step = threshold;

while ~feof(xml)
    
    if l > threshold; 
        
        disp(strcat('start loop line:',int2str(l))); 
        threshold = threshold + step;
    
    end
    
    line = fgetl(xml);
    l = l + 1;
    
end

fclose('all');

disp(strcat('Lines:',int2str(l)));

end