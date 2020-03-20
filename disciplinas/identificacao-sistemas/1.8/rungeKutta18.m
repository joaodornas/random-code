function rungeKutta = rungeKutta18(x,y,z,T,k)


xv = [];
yv = [];
zv = [];

xv(1) = x;
yv(1) = y;
zv(1) = z;

for i=2:k
    
    
    A = fxyz18( xv(i-1) , yv(i-1) , zv(i-1),'x');
    
    B = fxyz18( xv(i-1) + (T/2)*A , yv(i-1) , zv(i-1),'x');
   
    C = fxyz18( xv(i-1) + (T/2)*B , yv(i-1) , zv(i-1),'x');
    
    D = fxyz18( xv(i-1) + T*C , yv(i-1) , zv(i-1),'x');
    
    xv(i) = xv(i-1) + (T/6)*(A + 2*B + 2*C + D);
    
   
    A = fxyz18( xv(i-1) , yv(i-1) , zv(i-1),'y');
    
    B = fxyz18( xv(i-1) , yv(i-1) + (T/2)*A , zv(i-1),'y');
   
    C = fxyz18( xv(i-1) , yv(i-1)  + (T/2)*B , zv(i-1),'y');
    
    D = fxyz18( xv(i-1) , yv(i-1) + T*C , zv(i-1),'y');
    
    yv(i) = yv(i-1) + (T/6)*(A + 2*B + 2*C + D);
   
    
    A = fxyz18( xv(i-1) , yv(i-1) , zv(i-1),'z');
    
    B = fxyz18( xv(i-1) , yv(i-1) , zv(i-1) + (T/2)*A ,'z');
   
    C = fxyz18( xv(i-1) , yv(i-1) , zv(i-1) + (T/2)*B,'z');
    
    D = fxyz18( xv(i-1) , yv(i-1) , zv(i-1) + T*C ,'z');
    
    zv(i) = zv(i-1) + (T/6)*(A + 2*B + 2*C + D);
   
    
end

rungeKutta = [xv; yv; zv];
rungeKutta = rungeKutta.';


end