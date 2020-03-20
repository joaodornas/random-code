
prefix = 'Z:\Dropbox (Uni Magdeburg)\_DATA\LOW-HIGH-ATTENTION';
sufix = 'preprocessed\B0-DTI\6.forbedpost.bedpostX\Und_PRE_xx.mat';

DTI(1) = load(strcat(prefix,'\','SUBJECT-1-22-10-2015','\',sufix));
DTI(2) = load(strcat(prefix,'\','SUBJECT-2-26-10-2015','\',sufix));
DTI(3) = load(strcat(prefix,'\','SUBJECT-3-3-11-2015','\',sufix));
DTI(4) = load(strcat(prefix,'\','SUBJECT-4-2-11-2015','\',sufix));
DTI(5) = load(strcat(prefix,'\','SUBJECT-5-2-11-2015','\',sufix));
DTI(6) = load(strcat(prefix,'\','SUBJECT-6-24-11-2015','\',sufix));
DTI(7) = load(strcat(prefix,'\','SUBJECT-7-14-01-2016','\',sufix));
DTI(8) = load(strcat(prefix,'\','SUBJECT-8-14-01-2016','\',sufix));

aveC = zeros(size(DTI(1).C));

for iDTI=1:8
    
    aveC = aveC + DTI(iDTI).C;

end

aveC = aveC ./ 8;

DTI(9).C = aveC;

for iDTI=1:9
    
    for iiDTI=1:9
        
        cc(iDTI,iiDTI) = corr2(DTI(iDTI).C,DTI(iiDTI).C);
        
        if iDTI == iiDTI; cc(iDTI,iiDTI) = 0; end
        
    end
    
end

figure;

imagesc(cc);

min_C = 0.8;
max_C = 1;

caxis([min_C max_C]);
clrmp = colormap('jet');
clrmp(1,:) = [1 1 1];
colormap(clrmp);
colorbar;
%colorbar('Ticks',[min_C max_C/2 max_C]);

figure;

imagesc(aveC);

min_C = min(aveC(:));
max_C = max(aveC(:));

caxis([min_C max_C]);
clrmp = colormap('jet');
clrmp(1,:) = [1 1 1];
colormap(clrmp);
colorbar;
%colorbar('Ticks',[min_C max_C/2 max_C]);

