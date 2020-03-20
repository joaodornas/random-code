function h = movie2vector(videoin,video_name)


k = 1;

while (k < 301)

    framergb = read(videoin,k);
    
%     X = reshape(framergb.',[],1);
%     
%     X = X.';
     
    if k < 10
        frameindex = strcat('00', int2str(k));
    elseif k < 100
        frameindex = strcat('0', int2str(k));
    else
        frameindex = int2str(k);
    end

%     path = strcat(video_name,'_', frameindex);
% 
%     dlmwrite(strcat(path, '.txt'), X, 'delimiter', '\t', '-append', 'newline', 'pc');

%     imwrite(framergb,[video_name '-' frameindex '.jpg']);

    k = k + 1;

end
