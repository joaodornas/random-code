function hellinger_distance

Nsample = 1000;                                  % get two samples from same distribution

sample1 = gamrnd( 4, 10, 1, Nsample );

sample2 = gamrnd( 4, 10, 1, Nsample );



xbin = linspace(0, 1000, 501);                   % choose binning


p1 = hist( sample1, xbin ) / Nsample;            % probability vector of sample 1

p2 = hist( sample2, xbin ) / Nsample;            % probability vector of sample 2



HD = 1 - sum( sqrt( p1.*p2 ) )                   % compute Hellinger distance between original samples


sample1 = 5 * sample1;                           % change units / rescale values

sample2 = 5 * sample2;

p1 = hist( sample1, xbin ) / Nsample;   % probability vector of sample 1, rescaled

p2 = hist( sample2, xbin ) / Nsample;   % probability vector of sample 2, rescaled


HD = 1 - sum( sqrt( p1.*p2 ) )                   % compute Hellinger distance between rescaled samples


xbin = 5 * xbin;                                 % rescale binning
 


p1 = hist( sample1, xbin ) / Nsample;   % probability vector of sample 1, rescaled, with rescaled bins
p2 = hist( sample2, xbin ) / Nsample;   % probability vector of sample 2, rescaled, with rescaled bins


HD = 1 - sum( sqrt( p1.*p2 ) )                   % compute Hellinger distance between rescaled samples with rescaled bins


% conclusion: rescaling samples changes Hellinger distance if binning
% remains the same.   rescaling samples does not change Hellinger distance
% if binning is rescaled as well