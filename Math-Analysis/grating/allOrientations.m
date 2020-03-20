function allOrientations


how_many_directions = 1;

for i=0:22.5:360
    
    G = getGrating(200,0.95,125,2,0,i,57,how_many_directions);
    
    imwrite(G,strcat('grating-orientation-',num2str(i),'.jpg'),'jpeg');
    
    clear G;
    
end


end

