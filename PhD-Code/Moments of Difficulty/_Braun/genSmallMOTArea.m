

binsData = load('infoBinsCentersMOTArea.mat');

centers = binsData.info.centers;

prefixMOT = 'MOT-10-m-';
newprefixMOT = 'MOT-10-s-Area-';
prefixAnalyzeArea = 'MOT-Analyze-10-m-Area-';

nCenters = length(centers);
nMOT = length(centers(1).binsMOTidx);

sFrames = 3;   % sampling of frames
nWindow = 18/sFrames;  % number of sampled frames per window

FrameRate = 60;
Time = 10;

nTargets = 4;

for nC=1:nCenters
   
    for kMOT=1:nMOT
        
        filename = strcat(newprefixMOT,int2str(nC),'-',int2str(centers(nC).binsMOTidx(kMOT)),'.mat');
        
        originalMOT = load(strcat(prefixMOT,int2str(centers(nC).binsMOTidx(kMOT)),'.mat'));
        
        MOT = originalMOT.MOT;
        
        clear originalMOT
        
        kemin = centers(nC).binsMOTk(kMOT);
        
        kerange = kemin:kemin+nWindow;  % index to sampled frames;
        %jemin = (kemin-1)*sFrames+1;
        jerange = (kerange-1)*sFrames+1;  % index to original frames;

        start_ = jerange(1);
        end_ = jerange(end);
        
        start  = 61;               % skip first second
        finish = length(MOT.xi);   % 

        xi = MOT.xi(:,start:finish);
        yi = MOT.yi(:,start:finish);
        vi = MOT.vi(:,start:finish);
        vix = MOT.vix(:,start:finish);
        viy = MOT.viy(:,start:finish);
        ai = MOT.ai(:,start:finish);
        oi = MOT.oi(:,start:finish);
        
        framesBefore = Time*FrameRate/2;
        framesAfter = Time*FrameRate/2;
        
        start_ = start_ - framesBefore;
        end_ = end_ + framesAfter;
        
        range = start_:end_;
        
        MOT.xi = xi(:,range);
        MOT.yi = yi(:,range);
        MOT.vi = vi(:,range);
        MOT.vix = vix(:,range);
        MOT.viy = viy(:,range);
        MOT.ai = ai(:,range);
        MOT.oi = oi(:,range);
        
        MOT.AllRGB = MOT.Settings.MOT.Target.AllRGB;
                
        MOT.TargetColor = centers(nC).color(kMOT);
        
        MOT.Targets = centers(nC).targets(kMOT,:);
        
        MOT.Comb = centers(nC).comb(kMOT);
        
        MOT.MOTidx = centers(nC).binsMOTidx(kMOT);
        
        MOT.kidx = centers(nC).binsMOTk(kMOT);
        
        MOT.KArea = centers(nC).area(kMOT);

        analysis = load(strcat(prefixAnalyzeArea,int2str(centers(nC).binsMOTidx(kMOT)),'.mat'));
        
        frameBeforeInKs = round(start_/sFrames);
        frameAfterInKs = round(end_/sFrames);
        
        MOT.Areas = squeeze(analysis.areaComb(MOT.TargetColor,MOT.Comb,frameBeforeInKs:frameAfterInKs));
        
        clear analysis
        
        MOT.MomentOfDifficulty = 'Area';
        
        save(strcat(newprefixMOT,int2str(nC),'-',int2str(kMOT),'.mat'),'MOT');
        
    end
    
end
