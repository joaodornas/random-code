load('population_data','cells');
Sf_octave=-3:1;
Tf_octave=-2:4;

for i=11

save(strcat('/Users/lucaspinto/Documents/Lab/ProjectBooks/STTC-Book/Analyses/Multiple1DFits/workspaces/multiple_',char(cells(i))),'-append');

end