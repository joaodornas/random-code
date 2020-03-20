function calango(type)

start = 0;
finish = 10000;
start_time = 500;
end_time = 9500;

nFrames = 300;

limitTrial = 55;

pasta = '/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Calango/';
        
registros{1}.CR.x = 604;
registros{1}.CR.y = 402;
registros{1}.label = '_nsp006a01_1b-v3';
registros{1}.folder = '_nsp006a01_1b';
registros{1}.size.l = 44;
registros{1}.size.w = 37;
registros{1}.entrada.frame = 180;
registros{1}.entrada.time = registros{1}.entrada.frame * 30/1000 + start_time/1000;

registros{2}.CR.x = 598;
registros{2}.CR.y = 435;
registros{2}.label = '_nsp007a03_2a-v3';
registros{2}.folder = '_nsp007a03_2a';
registros{2}.size.l = 45;
registros{2}.size.w = 39;
registros{2}.entrada.frame = 215;
registros{2}.entrada.time = registros{2}.entrada.frame * 30/1000 + start_time/1000;

registros{3}.CR.x = 602;
registros{3}.CR.y = 362;
registros{3}.label = '_nsp008a02_1a-v3';
registros{3}.folder = '_nsp008a02_1a';
registros{3}.size.l = 25;
registros{3}.size.w = 18;
registros{3}.entrada.frame = 140;
registros{3}.entrada.time = registros{3}.entrada.frame * 30/1000 + start_time/1000;

registros{4}.CR.x = 569;
registros{4}.CR.y = 356;
registros{4}.label = '_nsp008a02_2b-v3';
registros{4}.folder = '_nsp008a02_2b';
registros{4}.size.l = 25;
registros{4}.size.w = 18;
registros{4}.entrada.frame = 132;
registros{4}.entrada.time = registros{4}.entrada.frame * 30/1000 + start_time/1000;

registros{5}.CR.x = 689;
registros{5}.CR.y = 442;
registros{5}.label = '_nsp011a01_1a-v3';
registros{5}.folder = '_nsp011a01_1a';
registros{5}.size.l = 42;
registros{5}.size.w = 30;
registros{5}.entrada.frame = 180;
registros{5}.entrada.time = registros{5}.entrada.frame * 30/1000 + start_time/1000;

registros{6}.CR.x = 683;
registros{6}.CR.y = 458;
registros{6}.label = '_nsp012a01_1b-v3';
registros{6}.folder = '_nsp012a01_1b';
registros{6}.size.l = 35;
registros{6}.size.w = 31;
registros{6}.entrada.frame = 184;
registros{6}.entrada.time = registros{6}.entrada.frame * 30/1000 + start_time/1000;

registros{7}.CR.x = 634;
registros{7}.CR.y = 271;
registros{7}.label = '_nsp033a09_1b-v311';
registros{7}.folder = '_nsp033a09_1b';
registros{7}.size.l = 42;
registros{7}.size.w = 28;
registros{7}.entrada.frame = 170;
registros{7}.entrada.time = registros{7}.entrada.frame * 30/1000 + start_time/1000;

registros{8}.CR.x = 634;
registros{8}.CR.y = 271;
registros{8}.label = '_nsp033a09_1c-v321';
registros{8}.folder = '_nsp033a09_1c';
registros{8}.size.l = 42;
registros{8}.size.w = 28;
registros{8}.entrada.frame = 170;
registros{8}.entrada.time = registros{8}.entrada.frame * 30/1000 + start_time/1000;

registros{9}.CR.x = 631;
registros{9}.CR.y = 273;
registros{9}.label = '_nsp033a09_3b-v331';
registros{9}.folder = '_nsp033a09_3b';
registros{9}.size.l = 42;
registros{9}.size.w = 28;
registros{9}.entrada.frame = 170;
registros{9}.entrada.time = registros{9}.entrada.frame * 30/1000 + start_time/1000;

Frame(1).start = 0;
Frame(1).end = 6;
Frame(1).label = 'parado';

Frame(2).start = 6;
Frame(2).end = 41;
Frame(2).label = 'andando';

Frame(3).start = 41;
Frame(3).end = 70;
Frame(3).label = 'parado';

Frame(4).start = 70;
Frame(4).end = 94;
Frame(4).label = 'andando';

Frame(5).start = 94;
Frame(5).end = 132;
Frame(5).label = 'parado';

Frame(6).start = 132;
Frame(6).end = 156;
Frame(6).label = 'andando';

Frame(7).start = 156;
Frame(7).end = 184;
Frame(7).label = 'andando';

Frame(8).start = 184;
Frame(8).end = 212;
Frame(8).label = 'andando';

Frame(9).start = 212;
Frame(9).end = 228;
Frame(9).label = 'parado';

Frame(10).start = 228;
Frame(10).end = 300;
Frame(10).label = 'andando';

%videoin = mmreader(['_MATRIZ-v' int2str(3) '-8bit-' int2str(1024) 'x' int2str(720) '.avi']);
    
%Frame = frameRF(Frame,videoin,registros,pasta);

%videoRF(registros,pasta);

for r=1:length(registros)
    
            %%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Load Spass...');

                Spass = load(strcat(char(registros{r}.label),'.mat'));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%% LOAD CONDITIONS TRIALS label   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Load Conditions Trials label...');

                nConditions = 1;
                
                for w=1:nConditions

                    trials_label(w).label = find(Spass.stimIds == w); 

                end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

            %%% SET RESOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Set Resolution...');

                spike_times = Spass.spike_times;

                spike_times = spike_times./ 32000;
                
                clear Spass;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('Calculate Spike Rate...');
            
            for w=1:nConditions

               trials_spikes = spike_times(trials_label(w).label,:);

               nTrials = size(trials_spikes,1);

               for q=1:nTrials

                   trial = trials_spikes(q,:);

                   trial = trial(trial>0);

                   trial = trial(trial>(start/1000) & trial<(finish/1000));
                   
                   for g=1:length(Frame)
                      
                       beginTime = ( Frame(g).start * 30/1000 ) + start_time/1000;
                       
                       %endTime = ( Frame(g).end * 30/1000 ) +
                       %start_time/1000;
                       
                       displacement = Frame(g).end - Frame(g).start;
                       
                       for l=1:displacement
                           
                            TrialRate(q).FrameRate(g).rate(l+1) = length(trial(trial>(beginTime + (l-1)*30/1000) & trial<(beginTime + l*30/1000))) / (30/1000);

                       end
                       
                       TrialRate(q).FrameRate(g).meanFrameRate = mean(TrialRate(q).FrameRate(g).rate);

                   end
           
                   registros{r}.condition(w).Trials(q).spikes = trial;
               
               end
               
                passo = limitTrial / 11;
                
                for g=1:length(Frame)
                    
                     displacement = Frame(g).end - Frame(g).start;
                    
                     for b=1:11

                            trialsCount = 0;

                            for x=((b-1)*passo + 1):(b*passo)

                                trial = trials_spikes(x,:);

                                trial = trial(trial>0);

                                trial = trial(trial>(start/1000) & trial<(finish/1000));

                                for l=1:displacement
                                                          
                                    trialsCount = trialsCount + length(trial(trial>(beginTime + (l-1)*30/1000) & trial<(beginTime + l*30/1000))) / (30/1000);

                                end
                                
                            end

                            registros{r}.frame(g).passo(b) = trialsCount / (displacement*passo);

                     end
                         
                end
                            
               
               for g=1:length(Frame)
                   
                   registros{r}.frame(g).rate = 0;
                   
                   for t=1:nTrials
                       
                        registros{r}.frame(g).rate = registros{r}.frame(g).rate + TrialRate(t).FrameRate(g).meanFrameRate;
                    
                   end
                   
                   registros{r}.frame(g).rate = registros{r}.frame(g).rate / nTrials;
                    
               end
                       
                range(1) = start/1000;
                range(2) = finish/1000;
                
                label = char(registros{r}.label);    
            
            end

end

save(strcat(pasta,'registros'),'registros');

for d=1:length(Frame)
   
    Frame(d).meanRate = 0;
    
    for r=1:length(registros)
        
        Frame(d).meanRate = Frame(d).meanRate + registros{r}.frame(d).rate; 
        
    end
    
    Frame(d).meanRate = Frame(d).meanRate / length(registros);
    
end

save(strcat(pasta,'Frame'),'Frame');

for r=1:length(registros)

    plot_raster(Frame,registros,r,range,label,start_time,type);        
    
end
                   

function plot_raster(Frame,registros,nRegistro,range,label,start_time,type)

    for k=1:length(registros{nRegistro}.condition)

      cla;
    
      hold on;
    
      idx=1;
        
      f = figure;
        
      for m = 1:length(registros{nRegistro}.condition(k).Trials)

        spikes = registros{nRegistro}.condition(k).Trials(m).spikes; 

        for i=1:length(spikes)

            plot([spikes(i) spikes(i)],[(idx-0.8) (idx)],'b');
            hold on;

        end

        idx = idx + 1;

      end

      for j=1:length(Frame)

          if mod(j,2) == 0
              
              color = 'r';
              
          else
              
              color = 'g';
              
          end
          
                p1 = [0 length(registros{nRegistro}.condition(k).Trials)];
                p2 = [(Frame(j).start * (30/1000) + start_time/1000) (Frame(j).start * (30/1000) + start_time/1000)];
                plot([p2(1) p2(1)],[p1(1) p1(2)],'Color',color,'LineWidth',0.5);
                hold on;
                
                if type == 1
                    
                    text(Frame(j).start*30/1000 + start_time/1000,30,int2str(registros{nRegistro}.frame(j).rate));
                    hold on;

                    text(Frame(j).start*30/1000 + start_time/1000,40,int2str(Frame(j).meanRate));
                    hold on;
                    
                elseif type == 2
                    
                    for o=1:11
                        
                        text(Frame(j).start*30/1000 + start_time/1000,5*o,int2str(registros{nRegistro}.frame(j).passo(o)));
                        hold on;
                        
                    end
                    
                end
                
      end
      
      p1 = [0 length(registros{nRegistro}.condition(k).Trials)];
      p2 = [registros{nRegistro}.entrada.time registros{nRegistro}.entrada.time];
      plot([p2(1) p2(1)],[p1(1) p1(2)],'Color','k','LineWidth',0.5);
      hold on;    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %hold off;
        %box on;
        xlabel('Tempo (s)');
        ylabel('{Repeti\c{c}\~oes}','interpreter','latex');
        %title(strrep(X.sites(n).label,'_','\_'));
        set(gca,'ylim',[0 idx-1]);
        set(gca,'ydir','rev');
        set(gca,'xlim',[range(1) range(2)]);
        set(gca,'xtick',unique([get(gca,'xtick') range]));

        folder = '/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Calango/';
        
        mkdir('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Calango',char(registros{nRegistro}.folder));
        
        path = strcat(folder,char(registros{nRegistro}.folder),'/');
        
        print(f,'-depsc',strcat(path,char(registros{nRegistro}.folder),'-','cond','-',int2str(k),'-','raster-walk-positions'));

        clear f;
        
    end
   
end

function Frame = frameRF(Frame,videoin,registros,pasta)
        
    for w=1:length(Frame)

        if w == 1,

            Frame(w).picture = read(videoin,1);

            Frame(w).x = start_time/1000;

        else

            Frame(w).picture = read(videoin,Frame(w).start);

            Frame(w).x = Frame(w).start * (30/1000) + start_time/1000;

        end

            for q=1:length(registros)

                    caminho = strcat(pasta,char(registros{q}.folder),'/');

                    X = registros{q}.CR.x;

                    Y = registros{q}.CR.y;

                    radius = round( max(registros{q}.size.l,registros{q}.size.w) / 2 );

                    yellow = uint8([255 255 0]);

                    circles = int32([X Y radius;X Y radius-1;X Y radius-2;X Y radius-3;X Y radius-4;X Y radius-5;X Y radius-6;X Y radius-7]); %;X Y radius-8;X Y radius-9;X Y radius-10;X Y radius-11;X Y radius-12;X Y radius-13;X Y radius-14]);

                    shapeInserter = vision.ShapeInserter('Shape','Circles','BorderColor','Custom','CustomBorderColor',yellow);

                    registros{q}.Frame(w).picture = step(shapeInserter,Frame(w).picture,circles);

                    mkdir('/Users/joaodornas/Documents/_Research/_EXPERIMENTOS/Forwards-Backwards/Experimento - v5/Calango',char(registros{q}.folder));
                    
                    imwrite(registros{q}.Frame(w).picture,strcat(caminho,char(registros{q}.folder),'-','Frame','-',int2str(w),'-','CR','.jpg'));

            end

    end
    
end

function videoRF(registros,pasta)
        
    for q=1:length(registros)

     caminho = strcat(pasta,char(registros{q}.folder),'/');

     videoObj(q) = VideoWriter(strcat(caminho,registros{q}.label,'-v3-RF-',int2str(registros{q}.size.l),'x',int2str(registros{q}.size.w),'.avi'));
     videoObj(q).FrameRate = 33;
     open(videoObj(q));

    end

    for q=1:length(registros)

            caminho = strcat(pasta,char(registros{q}.folder),'/');

            X = registros{q}.CR.x;

            Y = registros{q}.CR.y;

            radius = round( max(registros{q}.size.l,registros{q}.size.w) / 2 );

            yellow = uint8([255 255 0]);

            circles = int32([X Y radius;X Y radius-1;X Y radius-2;X Y radius-3;X Y radius-4;X Y radius-5;X Y radius-6;X Y radius-7]); %;X Y radius-8;X Y radius-9;X Y radius-10;X Y radius-11;X Y radius-12;X Y radius-13;X Y radius-14]);

            shapeInserter = vision.ShapeInserter('Shape','Circles','BorderColor','Custom','CustomBorderColor',yellow);

            for h=1:nFrames

                frame = read(videoin,h);

                frame = step(shapeInserter,frame,circles);

                writeVideo(videoObj(q),frame);

            end

    end


    for q=1:length(registros)

        close(videoObj(q));

    end
end

end

