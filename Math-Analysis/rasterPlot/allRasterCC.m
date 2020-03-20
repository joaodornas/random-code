function allRasterCC

registro{1} = '_nsp021a01_1a-v5.mat';
registro{2} = '_nsp021a02_1a-v6.mat';
registro{3} = '_nsp023a01_1a-v1.mat';
registro{4} = '_nsp023a02_1a-v4.mat';

for r=1:length(registro)

    rasterCC(registro{r})
    
    close all;
       
end

end

