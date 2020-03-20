function Grating = getGrating(dim,C,L,fs,ft,teta,distance,how_many_directions)


%SF = .25 .5 1 2 4 8 cycles/degree

%TF = .25 .5 1 2 4 8 Hz


minimum_view = 1;

visual_degree = round(2*atan(minimum_view/(2*distance)));


fs = fs * 0.0174532925;

teta = teta * (pi/180);
              
step_angle = 22.5 * (pi/180);

step_third = 30 * (pi/180);

anti = teta - pi;

semi_anti = teta - pi/2;

if anti < 0, anti = anti + 2*pi; end

if semi_anti < 0, semi_anti = semi_anti + 2*pi; end

switch how_many_directions
    
    case 1
        
        phi = [teta];
        
    case 2
        
        phi = [teta semi_anti];
            
    case 4
    
        phi = [teta (teta + step_third) (teta + 2*step_third) semi_anti];
    
    case 8
        
        phi = [teta (teta + step_angle) (teta + 2*step_angle) (teta + 3*step_angle) (teta + 4*step_angle) (teta + 5*step_angle) (teta + 6*step_angle) (teta + 7*step_angle) anti];
        
end

step_quadrante = 2*pi / how_many_directions;

for i=1:(length(phi))

    sumStepsQuadrantes(i) = i * step_quadrante;

end
   
diagonal = hypot(dim,dim);    
   
for x=1:dim

    for y=1:dim     
        
    switch how_many_directions
    
        case 1
        
            ypsilon = phi(1);
        
        case 2
        
            if (x > dim/2)
                
                ypsilon = phi(1);
                
            elseif (x <= dim/2)
                
                ypsilon = phi(2);
                
            end
                
            
        case 4
    
            if (x > dim/2)
                
                if (y > dim/2)
                    
                    ypsilon = phi(1);
                    
                elseif (y <= dim/2)
                    
                    ypsilon = phi(2);
                    
                end
                
            elseif (x <= dim/2)
                
                if (y <= dim/2)
                    
                    ypsilon = phi(3);
                    
                elseif (y > dim/2)
                    
                    ypsilon = phi(4);
                    
                end
                
            end
            
        end

        G(x,y) = L*( 1 + C*sin(2*pi*fs*(x*cos(ypsilon) + y*sin(ypsilon))-(2*pi*ft)) );

    end

 end
        
       
Grating = G;        

mask = imcircle(dim);

Grating((~mask)) = 30;

w = 4;
h = fspecial('gaussian',w,w);
Grating = imfilter(Grating,h);

Grating = uint8(Grating);



end

