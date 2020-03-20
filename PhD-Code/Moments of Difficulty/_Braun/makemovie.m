%makemovie(kECmin,kECmax,kEMmin,kEMmax,kEBmin,kEBmax,kEQmin,kEQmax,kVCmin,kVCmax,kVMmin,kVMmax,kVBmin,kVBmax,kVQmin,kVQmax,kVCballmin,kVCballmax,kVMballmin,kVMballmax,kVBballmin,kVBballmax,kVQballmin,kVQballmax,Ecolor,Emotion,Ebounce,Eshift,Lcolor,Lmotion,Lbounce,Lshift,Hcolor,Hmotion,Hbounce,Hshift,nFrames,sFrames,nWindow,nBalls,xi,yi,colorBalls,sizeX,sizeY,titlestring)

    close all;
    
    F = struct('cdata', [],'colormap', 'JET');
    
    videoObj = VideoWriter('entropies.avi');
    FrameRate = 60;
    open(videoObj);
    
    fs = 12;
    
    i = 1;
    
    for kemin=1:ceil(nFrames/sFrames)-nWindow;
        
        figure;
    
        set(gca,'nextplot','replacechildren');
    
        colormap jet;  
        
        kerange = kemin:kemin+nWindow-1;  % index to sampled frames;
        jemin = (kemin-1)*sFrames+1;
        jerange = (kerange-1)*sFrames+1;  % index to original frames;
        
        color{1} = 'r';
        color{2} = 'b';
        color{3} = 'k';

        for b=1:nBalls/2

            colorBalls(b) = 1;

        end
        for b=nBalls/2+1:nBalls

            colorBalls(b) = 2;
        end
        
        
        if kVCmin == kemin 
            
            plot( xi(kVCballmin,jerange(1)), yi(kVCballmin,jerange(1)), 'h', 'Color', color{colorBalls(kVCballmin)}, 'MarkerSize', 18 );
            colorBalls(kVCballmin) = 3; 
            
        end
        
        if kVCmax == kemin
            
            plot( xi(kVCballmax,jerange(1)), yi(kVCballmax,jerange(1)), 'h', 'Color', color{colorBalls(kVCballmax)}, 'MarkerSize', 18 );
            colorBalls(kVCballmax) = 3;
        
        end
        
        if kVMmin == kemin 
            
            plot( xi(kVMballmin,jerange(1)), yi(kVMballmin,jerange(1)), 'h', 'Color', color{colorBalls(kVMballmin)}, 'MarkerSize', 18 );
            colorBalls(kVMballmin) = 3;
        
        end
        
        if kVMmax == kemin 
            
            plot( xi(kVMballmax,jerange(1)), yi(kVMballmax,jerange(1)), 'h', 'Color', color{colorBalls(kVMballmax)}, 'MarkerSize', 18 );
            colorBalls(kVMballmax) = 3;
        
        end
        
        if kVBmin == kemin 
            
            plot( xi(kVBballmin,jerange(1)), yi(kVBballmin,jerange(1)), 'h', 'Color', color{colorBalls(kVBballmin)}, 'MarkerSize', 18 );
            colorBalls(kVBballmin) = 3; 
        
        end
        
        if kVBmax == kemin
            
            plot( xi(kVBballmax,jerange(1)), yi(kVBballmax,jerange(1)), 'h', 'Color', color{colorBalls(kVBballmax)}, 'MarkerSize', 18 );
            colorBalls(kVBballmax) = 3; 
        
        end
        
        if kVQmin == kemin 
            
            plot( xi(kVQballmin,jerange(1)), yi(kVQballmin,jerange(1)), 'h', 'Color', color{colorBalls(kVQballmin)}, 'MarkerSize', 18 );
            colorBalls(kVQballmin) = 3; 
        
        end
        
        if kVQmax == kemin 
            
            plot( xi(kVQballmax,jerange(1)), yi(kVQballmax,jerange(1)), 'h', 'Color', color{colorBalls(kVQballmax)}, 'MarkerSize', 18 );
            colorBalls(kVQballmax) = 3; 
        
        end
        
        
        hold on;
        for ib=1:nBalls

            plot( xi(ib,jerange), yi(ib,jerange), 'o', 'Color', color{colorBalls(ib)}, 'MarkerSize', 8 );
            plot( xi(ib,jerange(1)), yi(ib,jerange(1)), 'x', 'Color', color{3}, 'MarkerSize', 3 ); 

        end
        
        if kECmin == kemin; text( -50+(sizeX/2), 20, strcat('ECmin:',[num2str(Ecolor(kemin),'%4.2f')]), 'FontSize', fs ); end
        
        if kECmax == kemin; text( -50+(sizeX/2), 20, strcat('ECmax:',[num2str(Ecolor(kemin),'%4.2f')]), 'FontSize', fs ); end
        
        if kEMmin == kemin; text( -50+(sizeX/2), 40, strcat('EMmin:',[num2str(Emotion(kemin),'%4.2f')]), 'FontSize', fs ); end
        
        if kEMmax == kemin; text( -50+(sizeX/2), 40, strcat('EMmax:',[num2str(Emotion(kemin),'%4.2f')]), 'FontSize', fs ); end
        
        if kEBmin == kemin; text( -50+(sizeX/2), 60, strcat('EBmin:',[num2str(Ebounce(kemin),'%4.2f')]), 'FontSize', fs ); end
        
        if kEBmax == kemin; text( -50+(sizeX/2), 60, strcat('EBmax:',[num2str(Ebounce(kemin),'%4.2f')]), 'FontSize', fs ); end
        
        if kEQmin == kemin; text( -50+(sizeX/2), 80, strcat('EQmin:',[num2str(Eshift(kemin),'%4.2f')]), 'FontSize', fs ); end
        
        if kEQmax == kemin; text( -50+(sizeX/2), 80, strcat('EQmax:',[num2str(Eshift(kemin),'%4.2f')]), 'FontSize', fs ); end

        
        
        if kVCmin == kemin; text( -200+(sizeX/2), - 20, strcat('VCmin:',[num2str(Lcolor(kemin),'%4.2f')]), 'FontSize', fs ); end
        
        if kVCmax == kemin; text( -200+(sizeX/2), - 30, strcat('VCmax:',[num2str(Hcolor(kemin),'%4.2f')]), 'FontSize', fs ); end
        
        if kVMmin == kemin; text( -200+(sizeX/2), -20-20, strcat('VMmin:',[num2str(Lmotion(kemin),'%4.2f')]), 'FontSize', fs ); end
        
        if kVMmax == kemin; text( -200+(sizeX/2), -20-30, strcat('VMmax:',[num2str(Hmotion(kemin),'%4.2f')]), 'FontSize', fs ); end
        
        if kVBmin == kemin; text( -200+(sizeX/2), -40-20, strcat('VBmin:',[num2str(Lbounce(kemin),'%4.2f')]), 'FontSize', fs ); end
        
        if kVBmax == kemin; text( -200+(sizeX/2), -40-30, strcat('VBmax:',[num2str(Hbounce(kemin),'%4.2f')]), 'FontSize', fs ); end
        
        if kVQmin == kemin; text( -200+(sizeX/2), -60-20, strcat('VQmin:',[num2str(Lshift(kemin),'%4.2f')]), 'FontSize', fs ); end
        
        if kVQmax == kemin; text( -200+(sizeX/2), -60-30, strcat('VQmax:',[num2str(Hshift(kemin),'%4.2f')]), 'FontSize', fs ); end

     
        
        if kVCmin == kemin; text( -50+(sizeX/2), 0 - 20, strcat('ball:',int2str(kVCballmin)), 'FontSize', fs ); end
        
        if kVCmax == kemin; text( -50+(sizeX/2), 0 - 30, strcat('ball:',int2str(kVCballmax)), 'FontSize', fs ); end
        
        if kVMmin == kemin; text( -50+(sizeX/2), -20-20, strcat('ball:',int2str(kVMballmin)), 'FontSize', fs ); end
        
        if kVMmax == kemin; text( -50+(sizeX/2), -20-30, strcat('ball:',int2str(kVMballmax)), 'FontSize', fs ); end
        
        if kVBmin == kemin; text( -50+(sizeX/2), -40-20, strcat('ball:',int2str(kVBballmin)), 'FontSize', fs ); end
        
        if kVBmax == kemin; text( -50+(sizeX/2), -40-30, strcat('ball:',int2str(kVBballmax)), 'FontSize', fs ); end
        
        if kVQmin == kemin; text( -50+(sizeX/2), -60-20, strcat('ball:',int2str(kVQballmin)), 'FontSize', fs ); end
        
        if kVQmax == kemin; text( -50+(sizeX/2), -60-30, strcat('ball:',int2str(kVQballmax)), 'FontSize', fs ); end
             
        
        plot((-sizeX/2):(sizeX/2),0,'k');
        plot(0,(-sizeY/2):(sizeY/2),'k');
        
        text((-sizeX/2), -20-10, strcat('kemin:',int2str(kemin)), 'FontSize', fs);
        
        hold off;

        axis([-sizeX/2 sizeX/2 -sizeY/2 sizeY/2]);

        axis 'equal';
        axis 'off';
        
        F(kemin) = getframe;
 
        for i=1:FrameRate
            
            writeVideo(videoObj,F(kemin).cdata);  
            
        end

        close all;
        
        if mod((kemin*sFrames),60) == 0
            
            strcat(int2str(i),' second')
            i = i + 1;
            
        end
    
    end
    
    close(videoObj);
    
    save(strcat(titlestring,'.mat'),'F');
    
