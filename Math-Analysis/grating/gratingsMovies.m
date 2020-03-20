function gratingsMovies(dim,C,L,distance,Time_Exposition_in_Seconds,how_many_directions)

%SF = .25 .5 1 2 4 8 cycles/degree
%TF = .25 .5 1 2 4 8 Hz

FrameRate = 33.33;

%mainPath = strcat('/Volumes/Data/DATA/DTC/Grating/Videos/');

for ft = [0.25 0.5 1 2 4 8]

        for fs = [0.25 0.5 1 2 4 8]
            
             name = strcat('ft-',num2str(ft),'-fs-',num2str(fs),'-condicoes-16');
            
             writerObj = VideoWriter(strcat(name,'-video_out','.avi'));
             writerObj.FrameRate = FrameRate;
             open(writerObj);

            for teta=0:22.5:360
              
                Total_Frame_Rate = FrameRate*Time_Exposition_in_Seconds;
                
                ft_per_Frame = ft / FrameRate ; 
                
                ft_total = ft_per_Frame;
                
                for j=1:Total_Frame_Rate;

                    G = getGrating(dim,C,L,fs,ft_total,teta,distance,how_many_directions);
                    
                    writeVideo(writerObj,G);

                    ft_total = ft_total + ft_per_Frame;
                    
                end
                
            end
            
            close(writerObj);

        end
    
end

end

