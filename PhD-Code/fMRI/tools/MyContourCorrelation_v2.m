function MyContourCorrelation_v2(X,Y,label,compute_bins)

if compute_bins; getContourData(X,Y,label); end

plotContour(label);

end

function getContourData(X,Y,label)

% X is a vector of samples x_i, Y is a vector of samples y_i, both vectors have the same length


nsample = numel(X(:));

X = reshape( X, 1, nsample );

Y = reshape( Y, 1, nsample );

nbin = 20;

xbound = quantile( X, nbin-1 );

x_lo = [min(X) xbound];
x_hi = [xbound max(X)];
dx   = x_hi - x_lo;
x_mid = 0.5 * (x_lo + x_hi);

ybound = quantile( Y, nbin-1 );

y_lo = [min(Y) ybound];
y_hi = [ybound max(Y)];
dy   = y_hi - y_lo;
y_mid = 0.5 * (y_lo + y_hi);


D_xy = nan( nbin, nbin );  % joint density
P_xy = nan( nbin, nbin );  % joint probability
for i = 1 : nbin
    
    disp(strcat('bin:',int2str(i)));
    
    for j = 1 : nbin        
        idx = find( (x_lo(i) <= X) & (X <= x_hi(i)) & (y_lo(j) <= Y) & (Y <= y_hi(j)) );   % events in bin       
        P_xy( j, i )   = numel( idx ) / nsample; 
        D_xy( j, i )   = P_xy( j, i ) / ( dx(i) * dy(j) );
    end   
end

save(strcat('Contour-',label,'-bins.mat'),'x_mid','y_mid','D_xy','label');

end

function plotContour(label)

load(strcat('Contour-',label,'-bins.mat'));

f = figure;

hold on;
contour(x_mid, y_mid, D_xy);
contourcmap('jet','Colorbar','on');
plot([0 0.5],[0 0.5],'k-');
plot([0 1],[0 0.75],'k--');
if min(x_mid(:)) < 0; plot([0 -0.5],[0 0.5],'k-'); end
axis 'square';
xlabel('mean');
ylabel('STD');
if min(x_mid(:)) < 0; xlim([-1 1]); else xlim([0 1]); end
ylim([0 0.5]);
title(label);

print(f,'-depsc',strcat('Contour-',label,'.eps'));

end



