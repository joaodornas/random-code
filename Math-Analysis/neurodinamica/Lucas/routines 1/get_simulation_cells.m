cd '\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\PerroneWIM simulation'

%load workspaces containing names of elligible units

n_sust=size(sustCells,1);
n_trans=size(transCells,1);

%get sustained cell workspaces and write mean response and frequency
%matrices

for i=1:n_sust
    cellname=sustCells(i);
    filename=strcat('multiple_',char(cellname));
    load(filename,'STTCMeanResponse','Sf_octave','Tf_octave');
    assignin('base',strcat('sust_unit',int2str(i)),STTCMeanResponse')
    assignin('base',strcat('sust_unit',int2str(i),'_Sfs'),Sf_octave);
    assignin('base',strcat('sust_unit',int2str(i),'_Tfs'),Tf_octave);
end

%get transient cell workspaces and write mean response and frequency
%matrices

for i=1:n_trans
    cellname=transCells(i);
    filename=strcat('multiple_',char(cellname));
    load(filename,'STTCMeanResponse','Sf_octave','Tf_octave');
    assignin('base',strcat('trans_unit',int2str(i)),STTCMeanResponse')
    assignin('base',strcat('trans_unit',int2str(i),'_Sfs'),Sf_octave);
    assignin('base',strcat('trans_unit',int2str(i),'_Tfs'),Tf_octave);
end
   
save('\\.PSF\.Home\Documents\Lab\STTC-Book\Analyses\PerroneWIM simulation\simulation_cells')       