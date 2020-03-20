function MyContourCorrelation(X,Y,label)

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


figure;

% subplot(1,2,1);
% hold on;
% contour(x_mid, y_mid, P_xy);
% hold off;
% axis([x_lo(1) x_hi(end) y_lo(1) y_hi(end)]);
% axis 'square';
% xlabel('x');
% ylabel('y');
% title( 'Joint probability');

%subplot(1,2,2);
hold on;
contour(x_mid, y_mid, D_xy);
contourcmap('jet','Colorbar','on');
plot([0 0.5],[0 0.5],'k-');
plot([0 -0.5],[0 0.5],'k-');
%hold off;
% axis([x_lo(1) x_hi(end) y_lo(1) y_hi(end)]);
axis 'square';
xlabel('mean');
ylabel('STD');
xlim([-1 1]);
ylim([0 0.5]);
title(label);



