function h = viewSheilaData


data(1) = load('v3-12-08-20-sitio1-E1-res-15-bin-30-memoria-frame.mat');

data(2) = load('v3-12-08-27-sitio1-E1-res-15-bin-30-memoria-frame.mat');
data(3) = load('v3-12-08-27-sitio1-E2-res-15-bin-30-memoria-frame.mat');
data(4) = load('v3-12-08-29-sitio1-E1-res-15-bin-30-memoria-frame.mat');
data(5) = load('v4-12-08-29-sitio1-E1-res-15-bin-30-memoria-frame.mat');
data(6) = load('v1-12-08-29-sitio2-E1-res-15-bin-30-memoria-frame.mat');
data(7) = load('v1-12-08-30-sitio1-E1-res-15-bin-30-memoria-frame.mat');
data(8) = load('v4-12-08-30-sitio1-E1-res-15-bin-30-memoria-frame.mat');
data(9) = load('v1-12-08-30-sitio1-E3-res-15-bin-30-memoria-frame.mat');
data(10) = load('v3-12-08-30-sitio1-E3-res-15-bin-30-memoria-frame.mat');
data(11) = load('v4-12-08-30-sitio1-E3-res-15-bin-30-memoria-frame.mat');
data(12) = load('v1-12-08-31-sitio1-E1-res-15-bin-30-memoria-frame.mat');
data(13) = load('v3-12-08-31-sitio1-E1-res-15-bin-30-memoria-frame.mat');
data(14) = load('v4-12-08-31-sitio1-E1-res-15-bin-30-memoria-frame.mat');
data(15) = load('v1-12-08-31-sitio1-E3-res-15-bin-30-memoria-frame.mat');
data(16) = load('v3-12-08-31-sitio1-E3-res-15-bin-30-memoria-frame.mat');
data(17) = load('v4-12-08-31-sitio1-E3-res-15-bin-30-memoria-frame.mat');

data(18) = load('v2-12-09-03-sitio1-E3-res-15-bin-30-memoria-frame.mat');
data(19) = load('v4-12-09-03-sitio2-E1-res-15-bin-30-memoria-frame.mat');
data(20) = load('v4-12-09-04-sitio1-E1-res-15-bin-30-memoria-frame.mat');
data(21) = load('v3-12-09-05-sitio1-E1-res-15-bin-30-memoria-frame.mat');
data(22) = load('v4-12-09-05-sitio1-E1-res-15-bin-30-memoria-frame.mat');
data(23) = load('v3-12-09-06-sitio1-E1-res-15-bin-30-memoria-frame.mat');
data(24) = load('v4-12-09-06-sitio1-E1-res-15-bin-30-memoria-frame.mat');

    
   for L=1:size(data,2)
      
       maxL1(L) = data(L).memoryFrameData(1,1).maxL;
       maxL2(L) = data(L).memoryFrameData(1,2).maxL;
       maxL3(L) = data(L).memoryFrameData(1,3).maxL;
       
   end
   
maxL = struct('maxL1', maxL1, 'maxL2', maxL2, 'maxL3', maxL3);
    
save('maxL.mat','maxL');


end