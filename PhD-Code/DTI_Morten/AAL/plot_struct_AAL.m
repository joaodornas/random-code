
load('/Volumes/dropbox/_DATA/LOW-HIGH-ATTENTION/Jan-08-05-2015/preprocessed/B0-DTI/6.forbedpost.bedpostX/Und_PRE_xx.mat');

nROI = 90;

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 69 70];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [29 30 31 32 33 34 35 36 37 38 39 40 41 42 71 72 73 74 75 76 77 78];

x_idx(1) = 29;
x_idx(2) = 43;
x_idx(3) = 55;
x_idx(4) = 57;
x_idx(5) = 69;
x_idx(6) = 71;
x_idx(7) = 79;
x_idx(8) = 90;

for ix=1:length(x_idx)
    
    if ix == 1; x_idx_new(ix) = round(x_idx(ix)/2); end
    if ix ~= 1; x_idx_new(ix) = x_idx(ix-1) + round((x_idx(ix)-x_idx(ix-1))/2); end
    
end

x_label{1} = 'frontal';
x_label{2} = 'subcortical';
x_label{3} = 'occipital';
x_label{4} = 'temporal';
x_label{5} = 'parietal';
x_label{6} = 'frontal';
x_label{7} = 'subcortical';
x_label{8} = 'temporal';


% max_number_of_fibers = 5000;
% C = C./max_number_of_fibers;

imagesc(C);

max_C = max(C(:));
min_C = 0;

colorbar;
caxis([min_C max_C]);
colorbar('Ticks',[min_C max_C/2 max_C]);



for iROI=1:nROI-1
   
   if ~isempty(find(ismember(x_idx,iROI)))
       
       hold on;
       
       plot([iROI iROI],[1 nROI],'w-');
       
       plot([1 nROI],[iROI iROI],'w-');
       
   end
   
end

ax = gca;
ax.XTick = x_idx_new;
ax.XTickLabel = x_label;

ax.YTick = x_idx_new;
ax.YTickLabel = x_label;

xticklabel_rotate;