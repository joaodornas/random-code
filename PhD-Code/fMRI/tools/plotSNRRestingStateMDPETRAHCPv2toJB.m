function plotSNRRestingStateMDPETRAHCPv2

%%% SNR

% MD

load('Mean-STD-Voxels-AAL-MD.mat');

v(1).vector_one = std_voxels(:);
v(1).vector_two = mean_voxels(:);
v(1).title_label = 'MD';

clear mean_voxels std_voxels

% PETRA

load('Mean-STD-Voxels-AAL-PETRA.mat');

v(2).vector_one = std_voxels(:);
v(2).vector_two = mean_voxels(:);
v(2).title_label = 'PETRA';

clear mean_voxels std_voxels

% HCP

load('Mean-STD-Voxels-AAL-HCP.mat');

v(3).vector_one = std_voxels(:);
v(3).vector_two = mean_voxels(:);
v(3).title_label = 'HCP';

mysimpleplot(v);

clear mean_voxels std_voxels

end

function mysimpleplot(v)

color{1} = 'r*';
color{2} = 'b*';
color{3} = 'k*';

f = figure;

for iV=1:length(v)
    
    plot(v(iV).vector_two./v(iV).vector_one,v(iV).vector_one,color{iV},'MarkerSize',0.3);
    
    hold on;
    
end

legend({v(1).title_label,v(2).title_label,v(3).title_label});

end