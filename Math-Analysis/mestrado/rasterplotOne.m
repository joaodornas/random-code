function rasterplotOne(datapath,cur_color)

X = staread(strrep(datapath,'/',filesep));

  idx=1;
  for m=1:X.M
    for p = 1:X.categories(m).P
      for n = 1:X.N
        mins(idx) = X.categories(m).trials(p,n).start_time;
        maxes(idx) = X.categories(m).trials(p,n).end_time;
      end
    end
  end
  range(1) = min(mins);
  range(2) = max(maxes);
  range = X.sites(n).time_scale*range;

plot_raster(X,range,1,datapath,cur_color);



function plot_raster(X,range,n,datapath,cur_color)

f1 = figure;

idx=1;
cla;
hold on;
colororder = get(gca,'colororder');
num_colors = size(colororder,1);

for m=1:X.M

% switch m
%     
%     case 1
%         
%         cur_color = 'b';
%         
%     case 2
%         
%         cur_color = 'r';
%         
%     case 3
%         
%         cur_color = 'g';
%         
%     case 4
%         
%         cur_color = 'y';
%         
%     case 5
%         
%         cur_color = 'm';
%         
%     case 6
%         
%         cur_color = 'c';
%         
%     case 7
%         
%         cur_color = 'k';
%         
% end

  cur_start = idx;
  for p = 1:X.categories(m).P 
    list = X.sites(n).time_scale*X.categories(m).trials(p,n).list; 
    
    for i=1:size(list,2)
       
        plot([list(i) list(i)],[(idx-0.8) (idx)],cur_color);
        hold on;
    
    end
    
    idx = idx + 1;
  
  end
  
  cur_end = idx;

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

end

datapath = datapath(1:end-5);
print(f1,'-depsc',strcat(datapath,'-raster-plot-full'));

