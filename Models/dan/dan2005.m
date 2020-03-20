function dan2005

tic

disp('Load variables...');

registro = { 'nsp008a02_1a' };
channel = { 'E1' };
site_index = 1;
video_index = 3;
date = { '12-08-31' };

width = 1024;
height = 720;

x = 569;
y = 356;

radius = 68;

delay = 30;

start_time = 500;
end_time = 9500;
nFrames = 300;

nConditions = 2;

%%% LOAD VIDEO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Video...');

    videoin = mmreader(['_MATRIZ-v' int2str(video_index) '-8bit-' int2str(width) 'x' int2str(height) '.avi']);
    
    for i=1:nFrames
        
        frame = read(videoin,i);
        
        framergb(:,:) = frame(:,:,1);
       
        framecut = framergb((x-(radius/2)):(x+(radius/2)),(y-(radius/2)):(y+(radius/2)));
        
        movie(i,:,:) = framecut;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% LOAD SPASS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Spass...');

    Spass = load(strcat('_',char(registro{1}),'-','v',int2str(video_index),'.mat'));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%% LOAD CONDITIONS TRIALS LABELS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Load Conditions Trials Labels...');

    for i=1:nConditions
   
        trials_label(i).label = find(Spass.stimIds == i); 
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%% SET RESOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Set Resolution...');

    spike_times = Spass.spike_times;
    
    spike_times = spike_times./ 32000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% STC   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('STC...');

    CM = STC(spike_times,trials_label,nFrames,movie,delay,start_time,end_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%% EIGENVALUE AND EIGENVECTORS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Eigenvalues and Eigenvectors...');

    [V, D] = eig(CM);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%% WHITENING   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Whitening...');

    for i=1:nFrames
        
       WS(i,:,:) = movie(i,:,:).*V.*(((sqrt(D)).^-1));
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%%% STC   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('STC...');

    WC = STC(spike_times,trials_label,nFrames,WS,delay,start_time,end_time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%% EIGENVALUE AND EIGENVECTORS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Eigenvalues and Eigenvectors...');

    [WV, WD] = eig(WC);
    
    SV = (WV.').*(((sqrt(D)).^-1)).*(D)^1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


save('SV','SV');


function C = STC(spike_times,trials_label,nFrames,movie,delay,start_time,end_time)
%%% SPIKE-TRIGGERED COVARIANCE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Spike-Triggered Covariance...');

    forward = reshape(spike_times(trials_label(1).label,:).',[],1);
    
    forward = forward(forward>=start_time/1000 & forward<=end_time/1000);
    
    nSpikes = length(forward);
  
    for i=1:nFrames
    
        segment(:,:) = movie(i,:,:);
        
        S(i,:) = reshape(segment(:,:).',[],1);
        
    end
    
    C = zeros(size(S,2),size(S,2));
    
    for l=1:size(S,2)

        for w=1:size(S,2)

            for i=1:nSpikes

                timestamp = forward(i);
                
                stepbackms = timestamp - delay/1000;
                
                stepbackms = stepbackms - start_time/1000;
                
                stepback = round(stepbackms/(delay/1000));              
                
                if stepback > 0
                    
                    C(l,w) = C(l,w) + S(stepback,l)*S(stepback,w);
                    
                end

            end


        end

    end
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

end

        
end
    

