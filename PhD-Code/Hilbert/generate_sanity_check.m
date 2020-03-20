function sanity_check

%% create two concurrent sequences of time-points, one being the cause and the other being the effect

%  the purpose is to test the granger causality approach



Nevent  = 1000;  % exact number

Npoint  = 1000;  % AVERAGE number of time points per event


%  intervals between Poisson events

nuevent = 1/Npoint;

tisi    = -(1/nuevent) * log( rand(1,Nevent+2) );


%  event times and vector of time points

dummy  = cumsum(tisi);
Tevent = round(dummy(2:end-1));  % omit first and last time point

Tstart = 0;
Tend   = round(dummy(end));
Ti     = 1:Tend;


%  dominance values

Di     = -1*ones(size(Ti));

for i=1:2:Nevent
    
    idx_i = Tevent(i);
    
    if i<Nevent
        idx_f = Tevent(i+1)-1;
    else
        idx_f = Tend;
    end
    
    idx = idx_i:idx_f;
    
    Di(idx) = ones(size(idx));
    
end


%  cumulative history values

Ci     = zeros(size(Ti));

tauH  = 300;


for i=1:Nevent+1
    
    if i>1
        idx_i = Tevent(i-1);
    else
        idx_i = 1;
    end
    
    if i<=Nevent
        idx_f = Tevent(i);
    else
        idx_f = Tend;
    end
    
    idx = idx_i:idx_f;
    
    Cinf = Di(idx_i);
    
    C0   = Ci(idx_i);
    
    Ci(idx) = Cinf + ( C0 - Cinf ) * exp( -(idx - idx_i) / tauH );
       
end

size(Ti)
size(Di)
size(Ci)

% save values Ci and Di to file 

save( 'sanity.mat', 'Ci', 'Di', 'Tevent' );


% plot time sequences

figure;

hold on;
plot( Ti, Di, 'r.');
plot( Ti, Ci, 'b.');
hold off;

return;

