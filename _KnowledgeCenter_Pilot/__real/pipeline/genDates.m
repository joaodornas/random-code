clear all

start_year = 1900;
start_month = 01;
start_day = 01;

end_year = 2017;
end_month = 12;
end_day = 31;

first_day = datetime(start_year,start_month,start_day);
last_day = datetime(end_year,end_month,end_day);

period = first_day:last_day;

total_days = length(period);

matrix = zeros(total_days,3);

for iD=1:total_days
    
    this_day = period(iD);
    
    matrix(iD,1) = this_day.Day;
    matrix(iD,2) = this_day.Month;
    matrix(iD,3) = this_day.Year;
    
end

load('gcloud-zones.mat');

nZones = length(zones);
per_zone = 1000;

for iZ=1:nZones
   
    start_period = (1+(iZ-1)*per_zone);
    end_period = (per_zone+(iZ-1)*per_zone);
    
    if end_period > total_days
        
        end_period = total_days;
        
    end
    
    this_zone = zeros(length(start_period:end_period),3);
    
    this_zone(:,:,:) = matrix(start_period:end_period,:,:);
    
    dlmwrite(strcat(zones{iZ},'.txt'),this_zone,'delimiter',' ');
    
end
