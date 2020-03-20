
% load('all_primates_neuron.mat');
load('new_PMID.mat');
all_primates = new_PMID;
load('D:\__PDFs-Papers-for-Mac\RetinasPapersPerSource.mat');

count = 0;
for iPMID=1:length(EndNote)
    
    this_PMID = int2str(EndNote(iPMID));
    
    if ~isempty(find(strcmp(this_PMID,all_primates)))
       
        count = count + 1;
        gotEndNote{count} = this_PMID;
        
    end
    
end

count = 0;
for iPMID=1:length(ReadCube)
    
    this_PMID = int2str(ReadCube(iPMID));
    
    if ~isempty(find(strcmp(this_PMID,all_primates)))
       
        count = count + 1;
        gotReadCube{count} = this_PMID;
        
    end
    
end

count = 0;
for iPMID=1:length(PapersMac)
    
    this_PMID = int2str(PapersMac(iPMID));
    
    if ~isempty(find(strcmp(this_PMID,all_primates)))
       
        count = count + 1;
        gotPapersMac{count} = this_PMID;
        
    end
    
end

folder = 'D:\__PDFs-Papers-for-Mac\MISSING\PDF';
files = dir(strcat(folder,'\*.pdf'));

for iFile=1:length(files)

    MISSING{iFile} = files(iFile).name(1:end-4);

end

count = 0;
for iPMID=1:length(MISSING)
    
    this_PMID = MISSING(iPMID);
    
    if ~isempty(find(strcmp(this_PMID,all_primates)))
       
        count = count + 1;
        gotMISSING{count} = this_PMID;
        
    end
    
end

iigot = 0;
for igot=1:length(gotEndNote)
    iigot = iigot + 1;
    all_I_got(iigot) = str2num(gotEndNote{igot});
end
for igot=1:length(gotReadCube)
    iigot = iigot + 1;
    all_I_got(iigot) = str2num(gotReadCube{igot});
end
for igot=1:length(gotPapersMac)
    iigot = iigot + 1;
    all_I_got(iigot) = str2num(gotPapersMac{igot});
end
for igot=1:length(gotMISSING)
    iigot = iigot + 1;
    this_got = gotMISSING{igot};
    all_I_got(iigot) = str2num(this_got{1});
end

all_I_got = unique(all_I_got);

count = 0;
for iPMID=1:length(all_primates)

    this_PMID = all_primates(iPMID);

    if isempty(find(strcmp(this_PMID{1},all_I_got)))
       
        count = count + 1;
        stillMISSING{count} = this_PMID;
        
    end

end


