
function generateDocumentation(path)

mkdir(strcat(path,'/_doc'));

subDirectory = getDirs(path);
subDirectory = subDirectory';

iifun = 0;
for idir=1:size(subDirectory,1)
   
    subdir = subDirectory{idir,1};
    
    mfiles = dir(strcat(subdir,'/','*.m'));
    
    for imfile=1:length(mfiles)
        
       filename = mfiles(imfile).name;
       
       if ~strcmp(filename(1),'.')
       
            functions = fdep(filename);
            
            mlfun = functions.mlfun{1,1};
            
            if ~isempty(mlfun)
                
                for ifun=1:length(mlfun)
                    
                    iifun = iifun + 1;
                    
                    all_mlfun{iifun} = mlfun{ifun};
                
                end
                
            end
            
       end
        
    end
    
end

all_mlfun=all_mlfun';

[u_mlfun,c_mlfun] = count_unique(all_mlfun);

save(strcat(path,'/_doc','/','mfun.mat'),'u_mlfun','c_mlfun');

for iu=1:length(u_mlfun)
    
    funcfull = u_mlfun{iu,1};
    
    idx = strfind(funcfull,'/');
    
    funcname = funcfull(idx(end)+1:end);
    
    output{iu} = funcname;
    
end

fid = fopen(strcat(path,'/_doc/','mfun.txt'),'w');
for iout=1:length(output)
    fprintf(fid,'%s\n',output{iout});
end
fclose(fid);

end

function subDirectoryCellArr = getDirs(path)


    global subDirectoryCellArr
    global commandLineOutput
    global MAXNUMOFSUBDIRS
    
    MAXNUMOFSUBDIRS = 1000;
    
    directoryStr =    path;

    subDirectoryCellArr = cell(1,MAXNUMOFSUBDIRS);
    currCellIndex       = 0;
    commandLineOutput   = 1;
    
    % fill subDirectoryCellArr with all paths of subdirectories    
    currCellIndex = walkIn(directoryStr, currCellIndex);
    
end


function currCellIndex = walkIn(root, currCellIndex)    
    
    global subDirectoryCellArr
       
    if root(end) ~= filesep
        root(end+1) = filesep;
    end
    
    dirListing  = dir([root '*']);    
    
    for ii=1:length(dirListing)
        if dirListing(ii).isdir
            if dirListing(ii).name(1) ~= '.'
                % Add this folder path to listing of subdirectories
                currCellIndex                       = currCellIndex +1;
                subDirectoryCellArr{currCellIndex}  = [root dirListing(ii).name];
            
                % Look in this folder for more subdirectories
                currCellIndex = walkIn(subDirectoryCellArr{currCellIndex}, currCellIndex);                        
            end
        end
    end
end
