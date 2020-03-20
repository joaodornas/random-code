function k = fxyzRigid(x,y,z,s)


switch s
    
    case 'x'
   
        k = y * z;
        
    case 'y'
        
        k = -x * z;
        
    case 'z'
        
        k = -0.51 * x * y;
        
end


end