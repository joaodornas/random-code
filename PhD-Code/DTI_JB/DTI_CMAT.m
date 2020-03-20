
% prefix = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION';
% sufix = 'preprocessed\B0-DTI\6.forbedpost.bedpostX\Und_PRE_xx.mat';
% slashbar = '\';

prefix = '/Volumes/dropbox/_DATA/LOW-HIGH-ATTENTION';
sufix = 'preprocessed/B0-DTI/6.forbedpost.bedpostX/Und_PRE_xx.mat';
slashbar = '/';

nSubject = 8;

SUBJECT{1} = 'SUBJECT-1-22-10-2015';
SUBJECT{2} = 'SUBJECT-2-26-10-2015';
SUBJECT{3} = 'SUBJECT-3-3-11-2015';
SUBJECT{4} = 'SUBJECT-4-2-11-2015';
SUBJECT{5} = 'SUBJECT-5-2-11-2015';
SUBJECT{6} = 'SUBJECT-6-24-11-2015';
SUBJECT{7} = 'SUBJECT-7-14-01-2016';
SUBJECT{8} = 'SUBJECT-8-14-01-2016';

for iSubject=1:nSubject
    
    DTI(iSubject) = load([prefix slashbar SUBJECT{iSubject} slashbar sufix]);

end

for iSubject=1:nSubject
    
    dlmwrite(strcat('Connectome-',SUBJECT{iSubject},'.cmat'),DTI(iSubject).C,' ');
    
end