function VisualizeDTI

clear all;
close all;

%% load the matrices

Nsubject = 8;
Nsize    = 758;

MDTI = nan( Nsubject, Nsize, Nsize );

% for is = 1 : Nsubject
    for is = 1 : 5
        idx = is + 2;
    filename = ['DTI-Structural-Matrices/DTI-SUBJ-' num2str( idx, '%1d') '.mat']
    load( filename );
    MDTI(is,:,:) = C;
end

%% threshold the matrices

Cthreshold = 10;

Ntotal = sum(  MDTI(:) ~= 0 );

Ftotal = 100 * Ntotal / (Nsize * Nsize * Nsubject);

['fractional connectivity ' num2str( Ftotal, '%5.2f' ) ' %']

Nsubthresh = sum( MDTI(:) ~= 0 & MDTI(:) <= Cthreshold );

Fsubthresh = 100 * Nsubthresh / (Nsize * Nsize * Nsubject);

['subthreshold connectivity ' num2str( Fsubthresh, '%5.2f' ) ' %']

MDTI( MDTI ~= 0 & MDTI < Cthreshold ) = 0;

%% consistency of matrices

binarized_MDTI = ( MDTI ~= 0 );

summed_binarized_MDTI = squeeze( sum( binarized_MDTI, 1 ) );

vsum = 0:9;

[nsum, vsum] = hist( summed_binarized_MDTI(:), vsum );

nsum

% common connections

kcommon = find( summed_binarized_MDTI == 8 );

Fcommon = 100 * length( kcommon ) / (Nsize * Nsize);

fmedian_common = reshape( median( MDTI( :, kcommon ), 1 ), 1, length( kcommon(:) ) );

fvar_common    = reshape( var( MDTI( :, kcommon ), 1 ), 1, length( kcommon(:) ) );

[fmedian_common, idx] = sort( fmedian_common );
fvar_common    = fvar_common( idx );

% unique connections

kunique = find( summed_binarized_MDTI == 1 );

Funique = 100 * length( kunique ) / (Nsize * Nsize);

fmax_unique = reshape( max( MDTI( :, kunique )), 1, length( kunique(:) ) );

fvar_unique    = reshape( var( MDTI( :, kunique ), 1 ), 1, length( kunique(:) ) );

[fmax_unique, idx] = sort( fmax_unique );
fvar_unique    = fvar_unique( idx );


% partial connections

[fmedian_partial2, fvar_partial2] = PartialConnections( MDTI, 2 );
[fmedian_partial3, fvar_partial3] = PartialConnections( MDTI, 3 );
[fmedian_partial4, fvar_partial4] = PartialConnections( MDTI, 4 );
[fmedian_partial5, fvar_partial5] = PartialConnections( MDTI, 5 );
[fmedian_partial6, fvar_partial6] = PartialConnections( MDTI, 6 );
[fmedian_partial7, fvar_partial7] = PartialConnections( MDTI, 7 );

Ftwo   = 100 * length( fmedian_partial2 ) / (Nsize * Nsize);
Fthree = 100 * length( fmedian_partial3 ) / (Nsize * Nsize);
Ffour  = 100 * length( fmedian_partial4 ) / (Nsize * Nsize);
Ffive  = 100 * length( fmedian_partial5 ) / (Nsize * Nsize);
Fsix   = 100 * length( fmedian_partial6 ) / (Nsize * Nsize);
Fseven = 100 * length( fmedian_partial7 ) / (Nsize * Nsize);

% scatter plots median vs variance

figure;

subplot(2,2,1);
PlotMeanVar( fmedian_common, fvar_common, [1 0 0], 'common', fmax_unique, fvar_unique, [0 0 1], 'unique' );
title( [ num2str(Fcommon,'%5.2f') '% common / ' num2str(Funique,'%5.2f') '% unique'] );
subplot(2,2,2);
PlotMeanVar( fmedian_partial2, fvar_partial2, [0.4 0.4 0.4], 'two', fmedian_partial3, fvar_partial3, [0.2 0.2 0.8], 'three' );
title( [ num2str(Ftwo,'%5.2f') '% two / ' num2str(Fthree,'%5.2f') '% three'] );
subplot(2,2,3);
PlotMeanVar( fmedian_partial4, fvar_partial4, [0.4 0.4 0.4], 'four', fmedian_partial5, fvar_partial5, [0.2 0.8 0.2], 'five' );
title( [ num2str(Ffour,'%5.2f') '% four / ' num2str(Ffive,'%5.2f') '% five'] );
subplot(2,2,4);
hl = PlotMeanVar( fmedian_partial6, fvar_partial6, [0.4 0.4 0.4], 'six', fmedian_partial7, fvar_partial7, [0.8 0.2 0.2], 'seven' );
title( [ num2str(Fsix,'%5.2f') '% six / ' num2str(Fseven,'%5.2f') '% seven'] );


suptitle(  [ num2str(Ftotal,'%5.2f') '% total / ' num2str(Fsubthresh,'%5.2f') '% subthreshold'] );
axes(hl);

print 'DTIstats1' -depsc2;

% plot median vs density 

figure;

subplot(2,2,1);
PlotMeanDensity( fmedian_common, [1 0 0], 'common', fmax_unique, [0 0 1], 'unique' );
title( [ num2str(Fcommon,'%5.2f') '% common / ' num2str(Funique,'%5.2f') '% unique'] );
subplot(2,2,2);
PlotMeanDensity( fmedian_partial2, [0.4 0.4 0.4], 'two', fmedian_partial3, [0.2 0.2 0.8], 'three' );
title( [ num2str(Ftwo,'%5.2f') '% two / ' num2str(Fthree,'%5.2f') '% three'] );
subplot(2,2,3);
PlotMeanDensity( fmedian_partial4, [0.4 0.4 0.4], 'four', fmedian_partial5, [0.2 0.8 0.2], 'five' );
title( [ num2str(Ffour,'%5.2f') '% four / ' num2str(Ffive,'%5.2f') '% five'] );
subplot(2,2,4);
hl = PlotMeanDensity( fmedian_partial6, [0.4 0.4 0.4], 'six', fmedian_partial7, [0.8 0.2 0.2], 'seven' );
title( [ num2str(Fsix,'%5.2f') '% six / ' num2str(Fseven,'%5.2f') '% seven'] );


suptitle(  [ num2str(Ftotal,'%5.2f') '% total / ' num2str(Fsubthresh,'%5.2f') '% subthreshold'] );
axes(hl);

print 'DTIstats2' -depsc2;


%% output 7% and 16% connectivity matrices

binarized_MDTI = ( MDTI ~= 0 );

MDTI( ~binarized_MDTI ) = nan;   % replace zeros with nan

summed_binarized_MDTI = squeeze( sum( binarized_MDTI, 1 ) );

% 7% matrix with consensus connections
kk = find( summed_binarized_MDTI == 8 );
Mconsensus = nan( Nsize, Nsize );
Mconsensus(kk) = nanmedian( MDTI( :, kk ) );


% 16% matrix with majority connections
kk = find( summed_binarized_MDTI >= 5 );
Mmajority = nan( Nsize, Nsize );
Mmajority(kk) = nanmedian( MDTI( :, kk ) );

save( 'DTI_median_connectivity', 'Mconsensus', 'Mmajority' );

%% visualize two matrices

ShowConnectivity( Mconsensus );

title( 'Consensus connectivity, ca. 7%', 'FontSize', 14 );

print DTI_consensus_connectivity -dpdf;

ShowConnectivity( Mmajority );

title( 'Majority connectivity, ca. 16%', 'FontSize', 14 );

print DTI_majority_connectivity -dpdf;

return;

%% function defs

function [fmedian_partial, fvar_partial] = PartialConnections( MDTI, nconn )

[nobs, nsize, ~] = size( MDTI );

MDTI = reshape( MDTI, nobs, nsize*nsize );

size( MDTI )

binarized_MDTI = ( MDTI ~= 0 );

MDTI( ~binarized_MDTI ) = nan;   % replace zeros with nan

summed_binarized_MDTI = squeeze( sum( binarized_MDTI, 1 ) );

% kpartial = find( summed_binarized_MDTI == nconn );

fvalue_partial = MDTI( :, ( summed_binarized_MDTI == nconn ) );

fmedian_partial = nanmedian( fvalue_partial);    % compute mean over non-zero connections

fvar_partial = nanvar( fvalue_partial );

[fmedian_partial, idx] = sort( fmedian_partial );
fvar_partial    = fvar_partial( idx );

return;

%%

function hl = PlotMeanVar( fmean1, fvar1, clr1, lbl1, fmean2, fvar2, clr2, lbl2)

fs = 12;

hold on;
plot( -1, 0, 'o', 'Color', clr1);
plot( 0, -1, 'o', 'Color', clr2 );
hl = legend( lbl1, lbl2 );
set(hl,'Location','SouthEast', 'FontSize', fs );

plot( log10(fmean1), log10(sqrt(fvar1)), '.', 'LineWidth', 2, 'Color', clr1 );
plot( log10(fmean2), log10(sqrt(fvar2)), '.', 'LineWidth', 2, 'Color', clr2 );
hold off;

axis('square');
axis([1 4 0 3]);
xtic = (1:4);
xticlb = num2str(10.^(xtic'), '%1.0e');
ytic = (1:3);
yticlb = num2str(10.^(ytic'), '%1.0e');
set( gca, 'XTick', xtic, 'XTickLabel', xticlb, 'YTick', ytic, 'YTickLabel', yticlb, 'FontSize', fs );
xlabel('median','FontSize', fs);
ylabel('SD', 'FontSize', fs);


%%

function hl = PlotMeanDensity( fmean1, clr1, lbl1, fmean2, clr2, lbl2)

fs = 12;

hold on;
plot( -1, 0, 'o', 'Color', clr1);
plot( 0, -1, 'o', 'Color', clr2 );
hl = legend( lbl1, lbl2 );
set(hl,'Location','SouthEast', 'FontSize', fs );

n1 = length( fmean1 );
nbin = 101;
ix1 = round( linspace( 1, n1, nbin ) );
dens1 = zeros(size(ix1));
dens1(1) = (ix1(2)-ix1(1))/(n1*(fmean1(ix1(2))-fmean1(ix1(1))));
for i=2:length(ix1)-1
    dens1(i) = (ix1(i+1)-ix1(i-1))/(n1*(fmean1(ix1(i+1))-fmean1(ix1(i-1))));
end
dens1(end) = (ix1(end)-ix1(end-1))/(n1*(fmean1(ix1(end))-fmean1(ix1(end-1))));

n2 = length( fmean2 );
nbin = 101;
ix2 = round( linspace( 1, n2, nbin ) );
dens2 = zeros(size(ix1));
dens2(1) = (ix2(2)-ix2(1))/(n2*(fmean2(ix2(2))-fmean2(ix2(1))));
for i=2:length(ix2)-1
    dens2(i) = (ix2(i+1)-ix2(i-1))/(n2*(fmean2(ix2(i+1))-fmean2(ix2(i-1))));
end
dens2(end) = (ix2(end)-ix2(end-1))/(n2*(fmean2(ix2(end))-fmean2(ix2(end-1))));

plot( log10(fmean1(ix1)), dens1, '-', 'LineWidth', 2, 'Color', clr1 );
plot( log10(fmean2(ix2)), dens2, '-', 'LineWidth', 2, 'Color', clr2 );
hold off;

axis('square');
axis([1 4 0 0.1]);
xtic = (1:4);
xticlb = num2str(10.^(xtic'), '%1.0e');
ytic = (0.00:0.05:0.1);
yticlb = num2str(ytic', '%3.1f');
set( gca, 'XTick', xtic, 'XTickLabel', xticlb, 'YTick', ytic, 'YTickLabel', yticlb, 'FontSize', fs );
xlabel('median','FontSize', fs);
ylabel('density', 'FontSize', fs);

%%

function ShowConnectivity( Mconnectivity )


% prepate connectivity matrix
Mconnectivity = log( Mconnectivity );  % log scale of connection strength
cmin = min(Mconnectivity(:));
cmax = max(Mconnectivity(:));
Mconnectivity = round( 1.5 + 63 * ( Mconnectivity - cmin ) / (cmax - cmin) );  % scale [cmin cmax] to [2 64]
Mconnectivity( isnan( Mconnectivity ) ) = 1;  % assign 1 to empty connections

figure;

% prepare colormap
clrmp = colormap('jet');
clrmp(1,:) = [1 1 1];
clrmp(2:11,:) = repmat( [0.0 0.0 1.0], 10, 1);   % blue
clrmp(12:21,:) = repmat( [0.0 0.9 0.9], 10, 1);  % turqoise
clrmp(22:31,:) = repmat( [0 1.0 0.0], 10, 1);  % green
clrmp(32:41,:) = repmat( [1.0 1.0 0.0], 10, 1);  % yellow
clrmp(42:51,:) = repmat( [1.0 0.5 0.0], 10, 1);  % orange
clrmp(52:61,:) = repmat( [1.0 0.0 0.0], 10, 1);  % red
colormap( clrmp );
caxis( [0 cmax] );

% prepare anatomical labels
load('FC-Voxels-AAL-ROI-corr-KMeans-Info.mat');

ROI(23).label = 'Frontal-Sup-Med-L';
ROI(24).label = 'Frontal-Sup-Med-R';
ROI(69).label = 'Paracentral-Lob-L';
ROI(70).label = 'Paracentral-Lob-R';


nROI = 90;
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


xbox = 0.5 + [ROI_info{36,5} ROI_info{90,5} ROI_info{90,5} ROI_info{40,5} ROI_info{40,5} ROI_info{36,5} ROI_info{36,5}];
ybox = 0.5 + [ROI_info{36,5} ROI_info{36,5} ROI_info{40,5} ROI_info{40,5} ROI_info{90,5} ROI_info{90,5} ROI_info{36,5}];



hold on;

image( Mconnectivity );       % draw image

nClusters = ROI_info{nROI,5};

for i = 1 : length(xtic)      % draw anatomical boundaries
    
    jump = ROI_info{ xtic(i), 5};
    
    plot(0.5+jump*[1 1],0.5+[0 nClusters],'k-');
    plot([0.5+0 nClusters],0.5+jump*[1 1],'k-');
    
    if i>1
        xtic_jump(i) = 0.5 * ( ROI_info{ xtic(i-1), 5} + ROI_info{ xtic(i), 5} );
    else
        xtic_jump(i) = 0.5 * ( 1 + ROI_info{ xtic(i), 5} );
    end
end

plot( xbox, ybox, 'r-', 'LineWidth', 2 );

hold off;

ax = gca;
axis 'square';
axis([ 0.5 nClusters+0.5 0.5 nClusters+0.5] );
set(ax,'XTick', [] );
% set(ax,'XTick',xtic_jump);
% set(ax,'XTickLabel',xticlb);
%rotateXLabels( ax, 45 );

set(ax,'YTick',xtic_jump);
set(ax,'YTickLabel',xticlb);
%rotateYLabels( ax, 45 );

%xticklabel_rotate;

h = colorbar();             % make colorbar
ytic = [10 30 100 300 1000 3000];
yticlb = num2str( ytic' );
set( h, 'YTick', 1.5 + 63 * (log(ytic)-cmin) / (cmax-cmin), 'YTickLabel', yticlb );
ylabel( h, 'fibertracks (mean)', 'FontSize', 14 );

% print ROI labels on console

ROI_info
return;


