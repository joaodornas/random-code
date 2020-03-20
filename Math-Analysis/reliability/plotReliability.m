function plotReliability


registro = importdata('memoryBackwardProtocols.txt');

for r=1:length(registro)
       
    if r < 27
        
        protocol = registro{r}(2:end-7);
        
    else
        
        protocol = registro{r}(2:end-9);
        
    end
    
    for SIZE = 1:500
        
        getRealFor1D = load(strcat(protocol,'-reliabilityOptimal-','Forward','-','HSIZE','-',int2str(SIZE),'-','1D.mat'));
    
        %getRealBack1D = load(strcat(protocol,'-reliabilityOptimal-','Backward','-','HSIZE','-',int2str(SIZE),'-','1D.mat'));
        
        %getRealFor2D = load(strcat(protocol,'-reliabilityOptimal-','Forward','-','HSIZE','-',int2str(SIZE),'-','2D.mat'));
    
        %getRealBack2D = load(strcat(protocol,'-reliabilityOptimal-','Backward','-','HSIZE','-',int2str(SIZE),'-','2D.mat'));
    
    end
    
    f = figure;
    plot(1:500,getRealFor1D.inforFor.real_for,'r');
    hold on;
   
%     plot(1:500,getRealBack1D.inforBack.real_back,'b');
%     hold on;
%     
%     plot(1:500,getRealFor2D.inforFor.real_for,'y');
%     hold on;
%     
%     plot(1:500,getRealBack2D.inforBack.real_back,'k');
%     hold on;
%    
%     legend('Forward-1D','Backward-1D','Forward-2D','Backward-2D');
%     xlim([1 500]);
%     ylim([0 1]);
    
    text(400,0.9,strcat('real_for_max_1D:',num2str(max(getRealFor1D.inforFor.real_for))));
%     text(400,0.8,strcat('real_back_max_1D:',num2str(max(getRealBack1D.inforBack.real_back))));
%     text(400,0.7,strcat('real_for_max_2D:',num2str(max(getRealFor2D.inforFor.real_for))));
%     text(400,0.6,strcat('real_back_max_2D:',num2str(max(getRealBack2D.inforBack.real_back))));
    
    print(f,'-depsc',strcat('/Volumes/Data/DATA/Forward-Backward/reliabilityOptimal/',protocol,'-reliability-Backward-For&Back-1&2D.eps'));
   
    close all;
    
end




end
