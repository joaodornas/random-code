
clear all
clc

Pathway = '/Volumes/dropbox/_DATA/TVA/Pilot';

%Datfile = strcat(Pathway,'/CombiTVA_2346Targets-400-1.dat');
%Datfile = strcat(Pathway,'/CombiTVA_2346Targets-401-1.dat');
%Datfile = strcat(Pathway,'/CombiTVA_2346Targets-402-1.dat');
Datfile = strcat(Pathway,'/CombiTVA_2346Targets-403-1.dat');

tvadata = tvaloader(Datfile);

[theta,tvamodel,tvadata] = tvafit(tvadata);

tvareport(tvadata,tvamodel,theta);

[t,o,p,c] = tvaplot(tvadata,tvamodel,theta);

figure;
plot(t(1:6),o(1:6),'o');
hold on
plot(t(1:6),p(1:6),'-');



