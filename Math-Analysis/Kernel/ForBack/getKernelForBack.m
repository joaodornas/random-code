function kernel = getKernelForBack(registro,start_time,end_time,OS)

protocol = registro(2:end-4);

if strcmp(OS,'Mac')

    kernel = load(strcat('/Volumes/Data/DATA/Forward-Backward/kernel/',protocol,'-kernel-density-function-',int2str(start_time),'-',int2str(end_time),'.mat'));

elseif strcmp(OS,'Win')
    
   kernel = load(strcat('Z:\DATA\Forward-Backward\kernel\',protocol,'-kernel-density-function-',int2str(start_time),'-',int2str(end_time),'.mat'));
    
end

end

