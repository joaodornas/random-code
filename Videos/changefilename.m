

function h = changefilename(a)

for j=1:4
    
    for i=1:300
        
        if i < 10
            frame = strcat('00', int2str(i));
        elseif i < 100
            frame = strcat('0', int2str(i));
        else
            frame = int2str(i);
        end
        
        vj = strcat('/Users/joaodornas/Documents/videos-txt/v',int2str(j),'_backwards_720/','720x480-',frame,'.txt');
        vjd = strcat('/Users/joaodornas/Documents/videos-txt-d/v',int2str(j),'_backwards_720/v',int2str(j),'_backwards_720_',frame,'.txt');
        
        copyfile(vj,vjd);

    end
    
        for i=1:300
        
        if i < 10
            frame = strcat('00', int2str(i));
        elseif i < 100
            frame = strcat('0', int2str(i));
        else
            frame = int2str(i);
        end
        
        vj = strcat('/Users/joaodornas/Documents/videos-txt/v',int2str(j),'_backwards_1024/','1024x768-',frame,'.txt');
        vjd = strcat('/Users/joaodornas/Documents/videos-txt-d/v',int2str(j),'_backwards_1024/v',int2str(j),'_backwards_1024_',frame,'.txt');
        
        copyfile(vj,vjd);

        end
    
            for i=1:300
        
        if i < 10
            frame = strcat('00', int2str(i));
        elseif i < 100
            frame = strcat('0', int2str(i));
        else
            frame = int2str(i);
        end
        
        vj = strcat('/Users/joaodornas/Documents/videos-txt/v',int2str(j),'_forwards_720/','720x480-',frame,'.txt');
        vjd = strcat('/Users/joaodornas/Documents/videos-txt-d/v',int2str(j),'_forwards_720/v',int2str(j),'_forwards_720_',frame,'.txt');
        
        copyfile(vj,vjd);

            end
    
                for i=1:300
        
        if i < 10
            frame = strcat('00', int2str(i));
        elseif i < 100
            frame = strcat('0', int2str(i));
        else
            frame = int2str(i);
        end
        
        vj = strcat('/Users/joaodornas/Documents/videos-txt/v',int2str(j),'_forwards_1024/','1024x768-',frame,'.txt');
        vjd = strcat('/Users/joaodornas/Documents/videos-txt-d/v',int2str(j),'_forwards_1024/v',int2str(j),'_forwards_1024_',frame,'.txt');
        
        copyfile(vj,vjd);

    end
    
end

end

