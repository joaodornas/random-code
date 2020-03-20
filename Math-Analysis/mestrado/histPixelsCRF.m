function histPixelsCRF


pixelsCRF = [34, 27, 53, 30, 33, 31, 38, 32, 26, 52, 44, 35, 35, 48, 32, 50, 41, 43, 42, 35, 40, 38];
pixelsCRFgastc = [100, 68, 68, 136, 60, 68, 48, 40, 48, 68];

[upixelsCRF nupixelsCRF] = count_unique(pixelsCRF);

[upixelsCRFgastc nupixelsCRFgastc] = count_unique(pixelsCRFgastc);


[upixelsCRF sortIdx] = sort(upixelsCRF);
nupixelsCRF = nupixelsCRF(sortIdx);

[upixelsCRFgastc sortIdx] = sort(upixelsCRFgastc);
nupixelsCRFgastc = nupixelsCRFgastc(sortIdx);

frequenciasCRF = nupixelsCRF;
for i=2:length(frequenciasCRF)
   
    for k=1:i-1
    
        frequenciasCRF(i) = frequenciasCRF(i) + frequenciasCRF(k);
        
    end
        
end

frequenciasCRF = frequenciasCRF./max(frequenciasCRF);

frequenciasCRFgastc = nupixelsCRFgastc;
for i=2:length(frequenciasCRFgastc)
   
    for k=1:i-1
    
        frequenciasCRFgastc(i) = frequenciasCRFgastc(i) + frequenciasCRFgastc(k);
        
    end
        
end

frequenciasCRFgastc = frequenciasCRFgastc./max(frequenciasCRFgastc);

figure;

plot(upixelsCRF,frequenciasCRF,'b');
hold on;
plot(upixelsCRFgastc,frequenciasCRFgastc,'r');

% bar(histCRF);
% xlabel('Largura do CRC (pixels)','FontSize',30);
% ylabel('Quantidade Encontrada','FontSize',30);
% title('{Distribui\c{c}\~ao de Tamanhos de CRC GASTC}','interpreter','latex','FontSize',30);


% pixelsCRF([13 14 15 16 17 18 19 20 21 22]) = [];
% 
% median(pixelsCRF)
% median(pixelsCRFgastc)
% 
% h = lillietest(pixelsCRF)
% h = lillietest(pixelsCRFgastc)
% 
% [p,h] = ranksum(pixelsCRF,pixelsCRFgastc)

end

