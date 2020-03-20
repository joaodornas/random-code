
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

nTotalClusters = 758;
nClusters = 758;
nROI = 90;

hf = figure('color','w');
hold on;
hm = imagesc(ones(nClusters));

fs = 2;

% ROI and Clusters information
for n=1:nROI
    ROI_info{n,1} = n;                                        % AAL ROI number
    ROI_info{n,2} = ROI(n).label;                             % AAL ROI name
    ROI_info{n,3} = ROI(n).nClusters;                         % number of clusters in AAL ROI
    if n>1
        ROI_info{n,4} = ROI_info{n-1,5} + 1;                 % cumulative cluster number, first of range
        ROI_info{n,5} = ROI_info{n-1,5} + ROI(n).nClusters;  % cumulative cluster number, last of range
    else
        ROI_info{n,4} = 1;
        ROI_info{n,5} = ROI(n).nClusters;
    end
end

idx_frontal = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 69 70];
idx_rec_ins_cin = [27 28 29 30 31 32 33 34 35 36];
idx_hc_amyg = [37 38 39 40 41 42];
idx_occipital = [43 44 45 46 47 48 49 50 51 52 53 54];
idx_parietal = [57 58 59 60 61 62 63 64 65 66 67 68];
idx_temporal = [55 56 79 80 81 82 83 84 85 86 87 88 89 90];
idx_subcortical = [71 72 73 74 75 76 77 78];

cnt = 1;
xtic(cnt) = 2;
xticlb{cnt} = 'precent';
cnt=cnt+1;
xtic(cnt) = 6;
xticlb{cnt} = 'frontsup';
cnt=cnt+1;
xtic(cnt) = 10;
xticlb{cnt} = 'frontmid';
cnt=cnt+1;
xtic(cnt) = 16;
xticlb{cnt} = 'frontinf';
cnt=cnt+1;
xtic(cnt) = 18;
xticlb{5} = 'roland';
cnt=cnt+1;
xtic(cnt) = 20;
xticlb{cnt} = 'suppmot';
cnt=cnt+1;
xtic(cnt) = 22;
xticlb{cnt} = 'olfac';
cnt=cnt+1;
xtic(cnt) = 24;
xticlb{cnt} = 'frontsup';
cnt=cnt+1;
xtic(cnt) = 26;
xticlb{cnt} = 'frontmed';
cnt=cnt+1;
xtic(cnt) = 28;
xticlb{cnt} = 'rectus';
cnt=cnt+1;
xtic(cnt) = 30;
xticlb{cnt} = 'insula';
cnt=cnt+1;
xtic(cnt) = 36;
xticlb{cnt} = 'cingul';
cnt=cnt+1;
xtic(cnt) = 40;
xticlb{13} = 'HC';
cnt=cnt+1;
xtic(cnt) = 42;
xticlb{cnt} = 'amygd';
cnt=cnt+1;
xtic(cnt) = 44;
xticlb{cnt} = 'calcarine';
cnt=cnt+1;
xtic(cnt) = 46;
xticlb{cnt} = 'cuneus';
cnt=cnt+1;
xtic(cnt) = 48;
xticlb{cnt} = 'lingual';
cnt=cnt+1;
xtic(cnt) = 54;
xticlb{cnt} = 'occipit';
cnt=cnt+1;
xtic(cnt) = 56;
xticlb{cnt} = 'fusi';
cnt=cnt+1;
xtic(cnt) = 58;
xticlb{cnt} = 'postcent';
cnt=cnt+1;
xtic(cnt) = 62;
xticlb{cnt} = 'parietal';
cnt=cnt+1;
xtic(cnt) = 64;
xticlb{cnt} = 'suprmarg';
cnt=cnt+1;
xtic(cnt) = 66;
xticlb{cnt} = 'angul';
cnt=cnt+1;
xtic(cnt) = 68;
xticlb{cnt} = 'precun';
cnt=cnt+1;
xtic(cnt) = 70;
xticlb{cnt} = 'paracent';
cnt=cnt+1;
xtic(cnt) = 72;
xticlb{cnt} = 'caudate';
cnt=cnt+1;
xtic(cnt) = 74;
xticlb{cnt} = 'putamen';
cnt=cnt+1;
xtic(cnt) = 76;
xticlb{cnt} = 'pallidum';
cnt=cnt+1;
xtic(cnt) = 78;
xticlb{cnt} = 'thalamus';
cnt=cnt+1;
xtic(cnt) = 84;
xticlb{cnt} = 'tempsup';
cnt=cnt+1;
xtic(cnt) = 88;
xticlb{cnt} = 'tempmid';
cnt=cnt+1;
xtic(cnt) = 90;
xticlb{cnt} = 'tempinf';

face_alpha(1) = 0.1;
face_alpha(2) = 0.3;

face_color(1,:) = [125 125 125]./255;
face_color(2,:) = [200 200 200]./255;

for i = 1 : length(xtic)      % draw anatomical boundaries

    jump = ROI_info{ xtic(i), 5};

    if i == 1; 

        jump_before = 0; jump_after = jump; 

    else

        jump_before = ROI_info{ xtic(i-1), 5};
        jump_after = jump;

    end

%         plot(0.5+jump*[1 1],0.5+[0 nClusters],'k-');
%         plot([0.5+0 nClusters],0.5+jump*[1 1],'k-');

    for j=1:length(xtic)

         j_jump = ROI_info{ xtic(j), 5};

         if j == 1; 

            j_jump_before = 0; j_jump_after = j_jump; 

         else

            j_jump_before = ROI_info{ xtic(j-1), 5};
            j_jump_after = j_jump;

         end

         if mod(i,2) == 0

            if mod(j,2) == 0; a = 2; else a = 1; end

         else

             if mod(j,2) == 0; a = 1; else a = 2; end

         end

        fill([0.5+j_jump_before 0.5+j_jump_before 0.5+j_jump_after 0.5+j_jump_after],[0.5+jump_before 0.5+jump_after 0.5+jump_after 0.5+jump_before],face_color(a,:),'EdgeColor','none');

    end

    if i>1
        xtic_jump(i) = 0.5 * ( ROI_info{ xtic(i-1), 5} + ROI_info{ xtic(i), 5} );
    else
        xtic_jump(i) = 0.5 * ( 1 + ROI_info{ xtic(i), 5} );
    end
    
end

axis 'square';
axis([ 0.5 nClusters+0.5 0.5 nClusters+0.5] );

set(gca,'XTick',xtic_jump,'XTickLabel',xticlb,'YTick',xtic_jump,'YTickLabel',xticlb);
set(gca,'FontSize',fs);

xticklabel_rotate;

hold off;
 
print(hf,'-depsc','758-checkerboards.eps');
    
    