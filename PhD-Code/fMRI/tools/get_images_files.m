function images_files = get_images_files(fullpath,prefix_for_the_preprocessed_step,prefix,firstTR,lastTR)

images_files = char.empty;

iTR = 0;
for i=firstTR:lastTR
    
    iTR = iTR + 1;
    
    get_iTR_number;
    
    file = strcat(fullpath,'\',prefix_for_the_preprocessed_step,prefix,iTR_number,'-','0',iTR_number,'-','01','.hdr');
    
    images_files = strvcat(images_files,file);
            
end

end

