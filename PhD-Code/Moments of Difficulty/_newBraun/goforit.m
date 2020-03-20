function goforit

clear all;
close all;

v1 = -1+2*rand(2,200);

v2 = -1+2*rand(2,200);

figure;
Mkappa = pi/3;  
Msigma = pi/9; 
alpha = 0:0.01:pi;
efactor = exp( -(alpha-Mkappa) / Msigma );

Msimilarity = - (1 - efactor ) ./ ( 1 + efactor );
hold on;
plot( alpha*180/pi, Msimilarity, 'r.');
plot( 180*[alpha(1) alpha(end)]/pi, [0 0], 'k' );
plot( 180*[Mkappa Mkappa]/pi, [-1 1], 'k' );
hold off;



figure;

for i=1:20
    
    [A, B, C] = GetMsimilarity( v1(1,i), v1(2,i), v2(1,i), v2(2,i) );
    
    Malpha(i) = A;
    Mnorms(i) = B;
    Msimil(i) = C;
    
    hold off;
    plot( [0 0], [-1 1], 'k');
    axis 'equal';
    hold on;
    plot( [-1 1], [0 0], 'k');
    
    plot( [0 v1(1,i)], [0 v1(2,i)], 'r', 'LineWidth', 2 );
    plot( [0 v2(1,i)], [0 v2(2,i)], 'b', 'LineWidth', 2 );
    
    text( 0.5, 0.5, ['angle= ' num2str( Malpha(i)*180/pi, '%4.2f')], 'FontSize', 14 );
    text( 0.5, 0.4, ['norms= ' num2str( Mnorms(i), '%4.2f')], 'FontSize', 14 );
    text( 0.5, 0.3, ['simil= ' num2str( Msimil(i), '%4.2f')], 'FontSize', 14 );
    
    input('','s');


end
return;


function [Malpha, Mnorms, Msimilarity] = GetMsimilarity( vx1, vy1, vx2, vy2 )

    % measure motion similarity: 

    Mkappa = pi/3;  % threshold for motion similarity (smaller angles are similar, larger ones different )

    Msigma = pi/9;  % steepness of threshold

    Mroot = 0.25;  % choose exponent to match CV between alpha and norm distributions

    % get smaller angle between vectors

    alpha1 = atan2( vy1, vx1 );  % angle v1, from -pi to pi
    alpha2 = atan2( vy2, vx2 );  % angle v2, from -pi to pi
    alpha = alpha1 - alpha2;

    if alpha > pi
        alpha = alpha - 2*pi;
    end
    if alpha < -pi
        alpha = alpha + 2*pi;
    end
    Malpha = abs(alpha);
    
    efactor = exp( -(Malpha-Mkappa) / Msigma );

    Msimilarity = - (1 - efactor ) ./ ( 1 + efactor );

    [Malpha*180/pi (Malpha-Mkappa)*180/pi efactor Msimilarity]
    
    n1 = sqrt( vx1*vx1 + vy1*vy1 );  % norm 1

    n2 = sqrt( vx2*vx2 + vy2*vy2 );  % norm 2
    
    Mnorms = ( n1 * n2 )^Mroot;

    Msimilarity = Msimilarity * Mnorms;

return;