function stdMap = stdMap(stdM,y,x)

framestd = struct('cdata', zeros(y, x), 'colormap', gray);
framestd.cdata = stdM;

imwrite(framestd.cdata,'std-map.png');


end