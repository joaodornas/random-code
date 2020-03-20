function PICA_correlations_one_to_all

clear all;
close all;

global fs;

fs = 14;



load('Low-High-SUBJECT-1-groupICA-Four-lobes-All-Runs-AAL-residual-Parietal.mat');

% ICA_n(2)  number of PICA components (run1 and run2), per each of four lobes, and per each
% of 90 ROIs

% FC_Parietal(1).run(2).ROI(18) =
%
%     rho_high_occipital: [237x237 double]
%    pval_high_occipital: [237x237 double]
%     rho_rest_occipital: [237x237 double]
%    pval_rest_occipital: [237x237 double]
%      rho_high_parietal: [204x204 double]
%     pval_high_parietal: [204x204 double]
%      rho_rest_parietal: [204x204 double]
%     pval_rest_parietal: [204x204 double]
%      rho_high_temporal: [208x208 double]
%     pval_high_temporal: [208x208 double]
%      rho_rest_temporal: [208x208 double]
%     pval_rest_temporal: [208x208 double]
%       rho_high_frontal: [177x177 double]
%      pval_high_frontal: [177x177 double]
%       rho_rest_frontal: [177x177 double]
%      pval_rest_frontal: [177x177 double]
%

% FC_Parietal(1).ROI_Label =
%
% (1)  'Rolandic Oper L'
% (2)  'Rolandic Oper R'
% (3)  'Supp Motor Area L'
% (4)  'Supp Motor Area R'
% (5)  'Post Central L'
% (6)  'Post Central R'
% (7)  'Parietal Sup L'
% (8)  'Parietal Sup R'
% (9)  'Parietal Inf L'
% (10) 'Parietal Inf R'
% (11) 'Supramarginal L'
% (12) 'Supramarginal R'
% (13) 'Angular L'
% (14) 'Angular R'
% (15) 'Precuneus L'
% (16) 'Precuneus R'
% (17) 'Paracentral Lobule L'
% (18) 'Paracentral Lobule R'

Ncomponent = 40;

ROI_Number = 18;


for nPAIR = 9 : 2 : 17;
    
    nplot = 1;
    
    dLOBE     = zeros( 1, Ncomponent );
    dROI      = zeros( 1, Ncomponent );
    dLOBE_ROI = zeros( 1, Ncomponent );
    dROI_LOBE = zeros( 1, Ncomponent );
    
    ksLOBE     = zeros( 1, Ncomponent );
    ksROI      = zeros( 1, Ncomponent );
    ksLOBE_ROI = zeros( 1, Ncomponent );
    ksROI_LOBE = zeros( 1, Ncomponent );
    
    for nROI = nPAIR : nPAIR + 1;
        
        for nRUN = 1 : 2;
            
            %nROI = 18;
            
            
            ROI_ID = FC_Parietal(1).ROI_ID;
            
            ROIlabel = [char(FC_Parietal(1).ROI_Label( nROI )) ' Run' num2str(nRUN)];
            
            Nlobe = ICA_n( nRUN ).nICOccipital   % change lobe
            
            Nroi  = ICA_n( nRUN ).nICROI( ROI_ID( nROI ) )
            
            size( FC_Parietal(1).run( nRUN ).ROI( nROI ).rho_high_occipital )  % change lobe
            
            pairlabel = ['Occipital_ROI' num2str( ROI_ID( nROI ) ) ];  % change lobe
            
            rho_att = nan( 2*Ncomponent, 2*Ncomponent );
            pval_att = nan( 2*Ncomponent, 2*Ncomponent );
            rho_rest = nan( 2*Ncomponent, 2*Ncomponent );
            pval_rest = nan( 2*Ncomponent, 2*Ncomponent );
            
            for i=1:2
                ifrom = (i-1)*Nlobe      + [1:Ncomponent];
                ito   = (i-1)*Ncomponent + [1:Ncomponent];
                
                for j = 1:2
                    jfrom = (j-1)*Nlobe      + [1:Ncomponent];
                    jto   = (j-1)*Ncomponent + [1:Ncomponent];
                    
                    rho_att(ito,jto)   = FC_Parietal(1).run( nRUN ).ROI( nROI ).rho_high_occipital(ifrom,jfrom);   % change lobe
                    pval_att(ito,jto)  = FC_Parietal(1).run( nRUN ).ROI( nROI ).pval_high_occipital(ifrom,jfrom);
                    rho_rest(ito,jto)  = FC_Parietal(1).run( nRUN ).ROI( nROI ).rho_rest_occipital(ifrom,jfrom);
                    pval_rest(ito,jto) = FC_Parietal(1).run( nRUN ).ROI( nROI ).pval_rest_occipital(ifrom,jfrom);
                    
                end
            end
            
            ViewFC( ROIlabel, rho_att, pval_att, rho_rest, pval_rest, Ncomponent );
            
            %       print( [pairlabel '_' num2str(nplot) '_corrcoef'], '-depsc2' );
                    
            %        nplot = nplot + 1;
            
            [cLOBE, cROI, cLOBE_ROI, cROI_LOBE] = ViewNorm( ROIlabel, rho_att, pval_att, rho_rest, pval_rest, Ncomponent );
            
            %        print( [pairlabel '_norm'], '-depsc2' );
            
            %        ViewProjection( ROIlabel, rho_att, pval_att, rho_rest, pval_rest, Ncomponent )
            
            %        print( [pairlabel '_proj'], '-depsc2' );
            
            
            %dLOBE     = dLOBE     + NormContrast( cLOBE ); % cumulate CONTRAST between attention and resting
            %dROI      = dROI      + NormContrast( cROI );
            %dLOBE_ROI = dLOBE_ROI + NormContrast( cLOBE_ROI );
            %dROI_LOBE = dROI_LOBE + NormContrast( cROI_LOBE );
            
            dLOBE     = dLOBE     + NormDifference( cLOBE ); % cumulate DIFFERENCE between attention and resting
            dROI      = dROI      + NormDifference( cROI );
            dLOBE_ROI = dLOBE_ROI + NormDifference( cLOBE_ROI );
            dROI_LOBE = dROI_LOBE + NormDifference( cROI_LOBE );
            
            [kLOBE, kROI, kLOBE_ROI, kROI_LOBE] = DoKSTest( ROIlabel, rho_att, pval_att, rho_rest, pval_rest, Ncomponent );

            
            ksLOBE     = ksLOBE     + kLOBE;                 % cumulate REJECTIONS of KS null hypothesis
            ksROI      = ksROI      + kROI;
            ksLOBE_ROI = ksLOBE_ROI + kLOBE_ROI;
            ksROI_LOBE = ksROI_LOBE + kROI_LOBE;
            
        end
        

        
    end
    
    dLOBE     = dLOBE     / 4;
    dROI      = dROI      / 4;
    dLOBE_ROI = dLOBE_ROI / 4;
    dROI_LOBE = dROI_LOBE / 4;
    
    % view contrast
    
    % ViewNormContrast( ROIlabel(1:end-6), dLOBE, dROI, dLOBE_ROI, dROI_LOBE, Ncomponent );
    
    % print( [pairlabel '_norm_contrast'], '-depsc2' );
    
    % view difference
    
     ViewNormDifference( ROIlabel(1:end-6), dLOBE, dROI, dLOBE_ROI, dROI_LOBE, Ncomponent );
    
     print( [pairlabel '_norm_diff'], '-depsc2' );
     
     ViewKSResults( ROIlabel(1:end-6), ksLOBE, ksROI, ksLOBE_ROI, ksROI_LOBE, Ncomponent );
    
     print( [pairlabel '_KS_reject'], '-depsc2' );
   
    
    input('hit return to proceed');
        
    close all;
    
end

return;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nc = NormContrast( cMATRIX )

[nr, nc] = size( cMATRIX );

nc = zeros( 1, nc );

kk = find( cMATRIX(1,:) + cMATRIX(2,:) ~= 0 );  % avoid division by zero

nc(kk) = ( cMATRIX(1,kk) - cMATRIX(2,kk) ) ./ ( cMATRIX(1,kk) + cMATRIX(2,kk) );

return;

function nc = NormDifference( cMATRIX )

[nr, nc] = size( cMATRIX );

nc = zeros( 1, nc );

kk = find( cMATRIX(1,:) + cMATRIX(2,:) ~= 0 );  % avoid division by zero

nc(kk) = ( cMATRIX(1,kk) - cMATRIX(2,kk) );

return;

% csmatrix = zeros( 8, ROI_Number );
% 
% csmatrix( 1, : ) = csLOBE;
% csmatrix( 3, : ) = csLOBEROI;
% 
% csmatrix( 5, : ) = csROI;
% csmatrix( 7, : ) = csROILOBE;
% 
% sumlimit = 100;
% 
% summin = -sumlimit + sumlimit/32;
% summax =  sumlimit;
% sumrange = linspace( summin, summax, 64 );
% clrmp = colormap('jet');
% clrmp(32,:) = [1 1 1];
% 
% 
% 
% figure;
% 
% hold on;
% caxis([summin summax]);
% pcolor( csmatrix );
% colormap(clrmp);
% hold off;
% axis([0 ROI_Number+1 0 9]);
% 
% YTick = [1 3 5 7]+0.5;
% YTickLabel = { 'Lobe'; 'Lobe vs Roi'; 'Roi'; 'Roi vs Lobe'};
% 
% set( gca, 'YTick', YTick, 'YTickLabel', YTickLabel, 'FontSize', fs );
% 
% xlabel('ROI number', 'FontSize', fs );
% 
% 
% suptitle( titlestring );
% 
% print( titlestring, '-depsc2' );



return;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ViewFC( ROIlabel, rho_att, pval_att, rho_rest, pval_rest, Ncomponent )

global fs;

pcriterion = 0.001;
rholimit   = 0.5;

N2 = (length(rho_att(:)))^2;

% histogram correlations

figure;
subplot(2,2,1);
[N,X] = hist( rho_att(:) );
bar( X, N/N2 );
set(gca,'XLim',[-1 1]);
axis 'square';
xlabel( '\rho', 'FontSize', fs );
ylabel( 'attention', 'FontSize', fs );

subplot(2,2,2);
[N,X] = hist( pval_att(:) );
bar( X, N/N2 );
set(gca,'XLim',[-1 1]);
axis 'square';
xlabel( 'p-value', 'FontSize', fs );
ylabel( 'attention', 'FontSize', fs );

subplot(2,2,3);
[N,X] = hist( rho_rest(:) );
bar( X, N/N2 );
set(gca,'XLim',[-1 1]);
axis 'square';
xlabel( '\rho', 'FontSize', fs );
ylabel( 'resting', 'FontSize', fs );

subplot(2,2,4);
[N,X] = hist( pval_rest(:) );
bar( X, N/N2 );
set(gca,'XLim',[-1 1]);
axis 'square';
xlabel( 'p-value', 'FontSize', fs );
ylabel( 'resting', 'FontSize', fs );

suptitle(ROIlabel);

% prepare correlation plots

rhomin = -rholimit + rholimit/32;
rhomax =  rholimit;

rhorange = linspace( rhomin, rhomax, 64 );

kk = find( pval_att > pcriterion );
rho_att(kk) = zeros(size(kk));

kk = find( pval_rest > pcriterion );
rho_rest(kk) = zeros(size(kk));


clrmp = colormap('jet');
clrmp(32,:) = [1 1 1];

figure;
subplot(1,2,1);
hold on;
caxis([rhomin rhomax]);
h = pcolor( rho_att );
set( h, 'EdgeColor', 'none');
colormap(clrmp);
plot( 1+ Ncomponent*[1 1], Ncomponent*[0 2], 'k', 'LineWidth', 2 );
plot( Ncomponent*[0 2], 1+ Ncomponent*[1 1], 'k', 'LineWidth', 2 );
hold off;
axis 'equal'; 
Tick = [Ncomponent/2 3*Ncomponent/2];
TickLabel = { 'Lobe' ; 'ROI' };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );
title('attention', 'FontSize', fs );


h = colorbar;
set(h,'YLim',[-0.4 0.4], 'YTick',[-0.4 -0.2 0 0.2 0.4], 'PlotBoxAspectRatio', [1 20 1]);


subplot(1,2,2);
hold on;
caxis([rhomin rhomax]);
h = pcolor( rho_rest );
set( h, 'EdgeColor', 'none');
colormap(clrmp);
plot( 1+ Ncomponent*[1 1], Ncomponent*[0 2], 'k', 'LineWidth', 2 );
plot( Ncomponent*[0 2], 1+ Ncomponent*[1 1], 'k', 'LineWidth', 2 );
hold off;
axis 'equal'; 
Tick = [Ncomponent/2 3*Ncomponent/2];
TickLabel = { 'Lobe' ; 'ROI' };
set( gca, 'XTick', Tick, 'XTickLabel', TickLabel, 'YTick', Tick, 'YTickLabel', TickLabel );
title('resting', 'FontSize', fs );


h = colorbar;
set(h,'YLim',[-0.4 0.4], 'YTick',[-0.4 -0.2 0 0.2 0.4], 'PlotBoxAspectRatio', [1 20 1]);


suptitle( ROIlabel);

return;

%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cLOBE, cROI, cLOBE_ROI, cROI_LOBE] = ViewNorm( ROIlabel, rho_att, pval_att, rho_rest, pval_rest, Ncomponent )

global fs;

pcriterion = 0.001;
rholimit   = 0.5;


% prepare cumsum plots

rhomin = -rholimit + rholimit/32;
rhomax =  rholimit;

rhorange = linspace( rhomin, rhomax, 64 );

kk = find( pval_att > pcriterion );
rho_att(kk) = zeros(size(kk));

kk = find( pval_rest > pcriterion );
rho_rest(kk) = zeros(size(kk));


% prepare cumulative sums of significant corrcoefs

cROI  = nan(2, Ncomponent);  % ROI components, with each other, attention and resting
cLOBE = nan(2, Ncomponent);  % LOBE components, with each other, attention and resting

cROI_LOBE = nan(2, Ncomponent); % ROI components with LOBE components, attention and resting
cLOBE_ROI = nan(2, Ncomponent); % LOBE components with ROI components, attention and resting

for ic = 1 : Ncomponent
    
    % compute vector norms
    cLOBE(     1, ic )  = sqrt( nansum( abs( rho_att(               ic, 1           :  Ncomponent ) ).^2 ) );
    cLOBE(     2, ic )  = sqrt( nansum( abs( rho_rest(              ic, 1           :  Ncomponent ) ).^2 ) );
    cROI(      1, ic )  = sqrt( nansum( abs( rho_att(  Ncomponent + ic, Ncomponent+1:2*Ncomponent ) ).^2 ) );
    cROI(      2, ic )  = sqrt( nansum( abs( rho_rest( Ncomponent + ic, Ncomponent+1:2*Ncomponent ) ).^2 ) );
    
    cLOBE_ROI( 1, ic )  = sqrt( nansum( abs( rho_att(               ic, Ncomponent+1:2*Ncomponent ) ).^2 ) );
    cLOBE_ROI( 2, ic )  = sqrt( nansum( abs( rho_rest(              ic, Ncomponent+1:2*Ncomponent ) ).^2 ) );
    cROI_LOBE( 1, ic )  = sqrt( nansum( abs( rho_att(  Ncomponent + ic, 1           :  Ncomponent ) ).^2 ) );
    cROI_LOBE( 2, ic )  = sqrt( nansum( abs( rho_rest( Ncomponent + ic, 1           :  Ncomponent ) ).^2 ) );
    
    % compute sum of absolute components
    %cLOBE(     1, ic )  = nansum( abs( rho_att(               ic, 1           :  Ncomponent ) ) );
    %cLOBE(     2, ic )  = nansum( abs( rho_rest(              ic, 1           :  Ncomponent ) ) );
    %cROI(      1, ic )  = nansum( abs( rho_att(  Ncomponent + ic, Ncomponent+1:2*Ncomponent ) ) );
    %cROI(      2, ic )  = nansum( abs( rho_rest( Ncomponent + ic, Ncomponent+1:2*Ncomponent ) ) );
    
    %cLOBE_ROI( 1, ic )  = nansum( abs( rho_att(               ic, Ncomponent+1:2*Ncomponent ) ) );
    %cLOBE_ROI( 2, ic )  = nansum( abs( rho_rest(              ic, Ncomponent+1:2*Ncomponent ) ) );
    %cROI_LOBE( 1, ic )  = nansum( abs( rho_att(  Ncomponent + ic, 1           :  Ncomponent ) ) );
    %cROI_LOBE( 2, ic )  = nansum( abs( rho_rest( Ncomponent + ic, 1           :  Ncomponent ) ) );
    
end


figure;
subplot(1,2,1);
hold on;
plot( 1:Ncomponent, cLOBE(1,:), 'r:', 'LineWidth', 1 );
plot( 1:Ncomponent, cLOBE(2,:), 'b:', 'LineWidth', 1 );
plot( 1:Ncomponent, cLOBE_ROI(1,:), 'r', 'LineWidth', 2 );
plot( 1:Ncomponent, cLOBE_ROI(2,:), 'b', 'LineWidth', 2 );
hold off;
axis 'square';
xlabel( 'c LOBE', 'FontSize', fs );
ylabel( 'norm |cc|', 'FontSize', fs );

subplot(1,2,2);
hold on;
plot( 1:Ncomponent, cROI(1,:), 'r:', 'LineWidth', 1 );
plot( 1:Ncomponent, cROI(2,:), 'b:', 'LineWidth', 1 );
plot( 1:Ncomponent, cROI_LOBE(1,:), 'r', 'LineWidth', 2 );
plot( 1:Ncomponent, cROI_LOBE(2,:), 'b', 'LineWidth', 2 );
hold off;
axis 'square';
xlabel( 'c ROI', 'FontSize', fs );
ylabel( 'norm |cc|', 'FontSize', fs );


suptitle( ROIlabel);

return;

%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kLOBE, kROI, kLOBE_ROI, kROI_LOBE] = DoKSTest( ROIlabel, rho_att, pval_att, rho_rest, pval_rest, Ncomponent )

global fs;


pvalue = 0.05;       % pvalue for KS test



% prepare cumulative sums of significant corrcoefs

kROI  = nan(1, Ncomponent);  % ROI components, with each other, attention and resting
kLOBE = nan(1, Ncomponent);  % LOBE components, with each other, attention and resting

kROI_LOBE = nan(1, Ncomponent); % ROI components with LOBE components, attention and resting
kLOBE_ROI = nan(1, Ncomponent); % LOBE components with ROI components, attention and resting

for ic = 1 : Ncomponent
    
    % get samples
    x1            = rho_att(               ic, 1           :  Ncomponent );
    x2            = rho_rest(              ic, 1           :  Ncomponent );
    kLOBE(ic)     = kstest2( x1, x2, pvalue, 'larger' );
    
    x1            = rho_att(  Ncomponent + ic, Ncomponent+1:2*Ncomponent );
    x2            = rho_rest( Ncomponent + ic, Ncomponent+1:2*Ncomponent );
    kROI(ic)      = kstest2( x1, x2, pvalue, 'larger'  ); 
    
    x1            = rho_att(               ic, Ncomponent+1:2*Ncomponent );
    x2            = rho_rest(              ic, Ncomponent+1:2*Ncomponent );
    kLOBE_ROI(ic) = kstest2( x1, x2, pvalue, 'larger'  );
     
    x1            = rho_att(  Ncomponent + ic, 1           :  Ncomponent );
    x2            = rho_rest( Ncomponent + ic, 1           :  Ncomponent );
    kROI_LOBE(ic) = kstest2( x1, x2, pvalue, 'larger'  );
       
end

% compare cumulative distributions ...

N2 = Ncomponent*Ncomponent;

X1            = abs( rho_att(  1:Ncomponent, 1:Ncomponent ) );
X2            = abs( rho_rest( 1:Ncomponent, 1:Ncomponent ) );
CD1 = sort(X1(:));
CD2 = sort(X2(:));
ZD1 = sqrt(2) * erfinv( 2*CD1 - 1);
ZD2 = sqrt(2) * erfinv( 2*CD2 - 1);

figure;
subplot(2,2,1);
hold on;
plot( ZD1, [1:N2]/N2, 'r', 'LineWidth', 2 );
plot( ZD2, [1:N2]/N2, 'b', 'LineWidth', 2 );
hold off;
axis('square');
xlabel('z score');
ylabel('fraction');

title('LOBE-LOBE');

X1            = abs( rho_att(  Ncomponent+1:2*Ncomponent, Ncomponent+1:2*Ncomponent ) );
X2            = abs( rho_rest( Ncomponent+1:2*Ncomponent, Ncomponent+1:2*Ncomponent ) );
CD1 = sort(X1(:));
CD2 = sort(X2(:));
ZD1 = sqrt(2) * erfinv( 2*CD1 - 1);
ZD2 = sqrt(2) * erfinv( 2*CD2 - 1);

subplot(2,2,2);
hold on;
plot( ZD1, [1:N2]/N2, 'r', 'LineWidth', 2 );
plot( ZD2, [1:N2]/N2, 'b', 'LineWidth', 2 );
hold off;
axis('square');
xlabel('z score');
ylabel('fraction');

title('ROI-ROI');


X1            = abs( rho_att(  1:Ncomponent, Ncomponent+1:2*Ncomponent ) );
X2            = abs( rho_rest( 1:Ncomponent, Ncomponent+1:2*Ncomponent ) );
CD1 = sort(X1(:));
CD2 = sort(X2(:));
ZD1 = sqrt(2) * erfinv( 2*CD1 - 1);
ZD2 = sqrt(2) * erfinv( 2*CD2 - 1);

subplot(2,2,3);
hold on;
plot( ZD1, [1:N2]/N2, 'r', 'LineWidth', 2 );
plot( ZD2, [1:N2]/N2, 'b', 'LineWidth', 2 );
hold off;
axis('square');
xlabel('z score');
ylabel('fraction');

title('LOBE-ROI');

X1            = abs( rho_att(  Ncomponent+1:2*Ncomponent, 1:Ncomponent ) );
X2            = abs( rho_rest( Ncomponent+1:2*Ncomponent, 1:Ncomponent ) );
CD1 = sort(X1(:));
CD2 = sort(X2(:));
ZD1 = sqrt(2) * erfinv( 2*CD1 - 1);
ZD2 = sqrt(2) * erfinv( 2*CD2 - 1);

subplot(2,2,4);
hold on;
plot( ZD1, [1:N2]/N2, 'r', 'LineWidth', 2 );
plot( ZD2, [1:N2]/N2, 'b', 'LineWidth', 2 );
hold off;
axis('square');
xlabel('z score');
ylabel('fraction');

title('ROI-LOBE');

suptitle( ROIlabel );

return;

%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ViewProjection( ROIlabel, rho_att, pval_att, rho_rest, pval_rest, Ncomponent )

global fs;

pcriterion = 0.01;
rholimit   = 0.5;


% prepare cumsum plots

rhomin = -rholimit + rholimit/32;
rhomax =  rholimit;

rhorange = linspace( rhomin, rhomax, 64 );

kk = find( pval_att > pcriterion );
rho_att(kk) = zeros(size(kk));

kk = find( pval_rest > pcriterion );
rho_rest(kk) = zeros(size(kk));


% prepare projection of significant corrcoef

pROI  = nan(1, Ncomponent);  % ROI components with attention, projection on ROI components resting
pLOBE = nan(1, Ncomponent);  % LOBE components with attention, projection on LOBE components resting

pROI_LOBE = nan(1, Ncomponent); % ROI components with attention, projection on LOBE components resting
pLOBE_ROI = nan(1, Ncomponent); % LOBE components with attention, projection on ROI components resting

for ic = 1 : Ncomponent
    
    pLOBE(     1, ic )  = nansum( rho_att(               ic, 1           :  Ncomponent ) .* ...
                                  rho_rest(              ic, 1           :  Ncomponent ) );
    pROI(      1, ic )  = nansum( rho_att(  Ncomponent + ic, Ncomponent+1:2*Ncomponent ) .* ...
                                  rho_rest( Ncomponent + ic, Ncomponent+1:2*Ncomponent ) );
    
    pLOBE_ROI( 1, ic )  = nansum( rho_att(               ic, Ncomponent+1:2*Ncomponent ) .* ...
                                  rho_rest(              ic, Ncomponent+1:2*Ncomponent ) );
    pROI_LOBE( 1, ic )  = nansum( rho_att(  Ncomponent + ic, 1           :  Ncomponent ) .* ...
                                  rho_rest( Ncomponent + ic, 1           :  Ncomponent ) );
    
end


figure;
subplot(1,2,1);
hold on;
plot( 1:Ncomponent, pLOBE(1,:), 'b--', 'LineWidth', 2 );
plot( 1:Ncomponent, pLOBE_ROI(1,:), 'r', 'LineWidth', 2 );
h = legend( 'on LOBE', 'on ROI');
set( h, 'FontSize', fs);
plot( [0 Ncomponent], [0 0], 'k', 'LineWidth', 1 );
hold off;
axis 'square';
xlabel( 'c LOBE', 'FontSize', fs );
ylabel( 'proj cc', 'FontSize', fs );

subplot(1,2,2);
hold on;
plot( 1:Ncomponent, pROI(1,:), 'b--', 'LineWidth', 2 );
plot( 1:Ncomponent, pROI_LOBE(1,:), 'r', 'LineWidth', 2 );
h = legend( 'on ROI', 'on LOBE');
set( h, 'FontSize', fs);
plot( [0 Ncomponent], [0 0], 'k', 'LineWidth', 1 );
hold off;
axis 'square';
xlabel( 'c ROI', 'FontSize', fs );
ylabel( 'proj cc', 'FontSize', fs );


suptitle( ROIlabel);

return;

%%

function ViewNormContrast( ROIlabel, dLOBE, dROI, dLOBE_ROI, dROI_LOBE, Ncomponent )

global fs;

hmax = 0.5;

figure;
subplot(1,2,1);
hold on;
plot( 1:Ncomponent, dLOBE, 'b:', 'LineWidth', 2 );
plot( 1:Ncomponent, dLOBE_ROI, 'r', 'LineWidth', 2 );
h = legend( 'with LOBE cs', 'with ROI cs');
set( h, 'FontSize', fs );
plot( [0 Ncomponent], [1 1]*nanmean(dLOBE), 'b', 'LineWidth', 1 );
plot( [0 Ncomponent], [1 1]*nanmean(dLOBE_ROI), 'r', 'LineWidth', 1 );
plot( [0 Ncomponent], [0 0], 'k', 'LineWidth', 1 );
hold off;
axis 'square';
axis([0 Ncomponent -hmax hmax]);
xlabel( 'LOBE c', 'FontSize', fs );
ylabel( 'ave diff', 'FontSize', fs );

subplot(1,2,2);
hold on;
plot( 1:Ncomponent, dROI, 'b:', 'LineWidth', 2 );
plot( 1:Ncomponent, dROI_LOBE, 'r', 'LineWidth', 2 );
h = legend( 'with ROI cs', 'with LOBE cs');
set( h, 'FontSize', fs );
plot( [0 Ncomponent], [1 1]*nanmean(dROI), 'b', 'LineWidth', 1 );
plot( [0 Ncomponent], [1 1]*nanmean(dROI_LOBE), 'r', 'LineWidth', 1 );
plot( [0 Ncomponent], [0 0], 'k', 'LineWidth', 1 );
hold off;
axis 'square';
axis([0 Ncomponent -hmax hmax]);
xlabel( 'd ROI', 'FontSize', fs );
ylabel( 'ave diff', 'FontSize', fs );


suptitle( ROIlabel);

return;

%%

function ViewNormDifference( ROIlabel, cLOBE, cROI, cLOBE_ROI, cROI_LOBE, Ncomponent )

global fs;

hmax = 0.5;

figure;
subplot(1,2,1);
hold on;
plot( 1:Ncomponent, cLOBE, 'b:', 'LineWidth', 2 );
plot( 1:Ncomponent, cLOBE_ROI, 'r', 'LineWidth', 2 );
h = legend( 'with LOBE cs', 'with ROI cs');
set( h, 'FontSize', fs );
plot( [0 Ncomponent], [1 1]*nanmean(cLOBE), 'b', 'LineWidth', 1 );
plot( [0 Ncomponent], [1 1]*nanmean(cLOBE_ROI), 'r', 'LineWidth', 1 );
plot( [0 Ncomponent], [0 0], 'k', 'LineWidth', 1 );
hold off;
axis 'square';
axis([0 Ncomponent -hmax hmax]);
xlabel( 'LOBE c', 'FontSize', fs );
ylabel( 'ave diff', 'FontSize', fs );

subplot(1,2,2);
hold on;
plot( 1:Ncomponent, cROI, 'b:', 'LineWidth', 2 );
plot( 1:Ncomponent, cROI_LOBE, 'r', 'LineWidth', 2 );
h = legend( 'with ROI cs', 'with LOBE cs');
set( h, 'FontSize', fs );
plot( [0 Ncomponent], [1 1]*nanmean(cROI), 'b', 'LineWidth', 1 );
plot( [0 Ncomponent], [1 1]*nanmean(cROI_LOBE), 'r', 'LineWidth', 1 );
plot( [0 Ncomponent], [0 0], 'k', 'LineWidth', 1 );
hold off;
axis 'square';
axis([0 Ncomponent -hmax hmax]);
xlabel( 'd ROI', 'FontSize', fs );
ylabel( 'ave diff', 'FontSize', fs );


suptitle( ROIlabel);

return;

%%

function ViewKSResults( ROIlabel, ksLOBE, ksROI, ksLOBE_ROI, ksROI_LOBE, Ncomponent )

global fs;

hmax = 4;

figure;
subplot(1,2,1);
hold on;
stairs( 1:Ncomponent, ksLOBE, 'b:', 'LineWidth', 2 );
stairs( 1:Ncomponent, ksLOBE_ROI, 'r', 'LineWidth', 2 );
h = legend( 'with LOBE cs', 'with ROI cs');
set( h, 'FontSize', fs );
hold off;
axis 'square';
axis([0 Ncomponent -hmax hmax]);
xlabel( 'LOBE c', 'FontSize', fs );
ylabel( 'KS reject', 'FontSize', fs );

subplot(1,2,2);
hold on;
stairs( 1:Ncomponent, ksROI, 'b:', 'LineWidth', 2 );
stairs( 1:Ncomponent, ksROI_LOBE, 'r', 'LineWidth', 2 );
h = legend( 'with ROI cs', 'with LOBE cs');
set( h, 'FontSize', fs );
hold off;
axis 'square';
axis([0 Ncomponent -hmax hmax]);
xlabel( 'ROI c', 'FontSize', fs );
ylabel( 'KS reject', 'FontSize', fs );


suptitle( ROIlabel);

return;


