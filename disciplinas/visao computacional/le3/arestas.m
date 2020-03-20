function outimgrgb = arestas(inimgpath, N)

img = imread(inimgpath);

[L1, L2, V1] = corners(img, N);

outimg = L1;

il = size(img,1);
jl = size(img,2);

angles = [0; 45; 90; 135];

maxvalueL1 = max(max(L1));
maxvalueL2 = max(max(L2));

for i=1:il
    
   for j=1:jl
   
       orientation = abs(atand(V1(i,j,2)/V1(i,j,1)));
       [minDifAngle, indexAngle] = min(abs(angles-orientation));
       
       
       switch indexAngle
           
           case 1
      
                ia = i;
                id = i;
                ja = j - 1;
                jd = j + 1;
      
           case 2
      
                ia = i + 1;
                id = i - 1;
                ja = j - 1;
                jd = j + 1;
      
           case 3
      
                ia = i - 1;
                id = i + 1;
                ja = j;
                jd = j;
      
           case 4
      
                ia = i - 1;
                id = i + 1;
                ja = j - 1;
                jd = j + 1;
      
       end
          
           
       if ia < 1
           ia = 1;
       elseif ia > il
           ia = il;
       end
       if id < 1
           id = 1;
       elseif id > il
           id = il;
       end
       if ja < 1
           ja = 1;
       elseif ja > jl
           ja = jl;
       end
       if jd < 1
           jd = 1;
       elseif jd > jl
           jd = jl;
       end
       
    
       if ((L1(i,j) > L1(ia,ja)) && (L1(i,j) > L1(id,jd)) && (L1(i,j) > (0.005*maxvalueL1)))
          
           outimg(i,j) = 0;
           
       else
           
           outimg(i,j) = 255;
           
       end
     
   end
    
end

outimgrgb = cat(3,outimg,outimg,outimg);

for i=1:il
    
    for j=1:jl
       
       if (L1(i,j) >= L2(i,j) && L2(i,j) > (0.85*maxvalueL2) && outimg(i,j) == 0)
             
           N = 10;
           
           for iq=0:(2*N)
              
              for jq=0:(2*N)
                 
                  if (iq == 0 || iq == 2*N)
                      
                      outimgrgb(i+iq,j+jq,1) = 255;

                      outimgrgb(i+iq,j+jq,2) = 0;

                      outimgrgb(i+iq,j+jq,3) = 0;
                  
                  elseif (jq == 0 || jq == 2*N)  
                      
                      outimgrgb(i+iq,j+jq,1) = 255;

                      outimgrgb(i+iq,j+jq,2) = 0;

                      outimgrgb(i+iq,j+jq,3) = 0;
                      
                  end
                  
              end
              
           end
        
       end
       
    end
    
end

end

