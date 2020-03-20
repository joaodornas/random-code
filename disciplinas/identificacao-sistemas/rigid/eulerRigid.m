function euler = eulerRigid(x,y,z,T,k)

xv = [];
yv = [];
zv = [];

xv(1) = x;
yv(1) = y;
zv(1) = z;

for i=2:k
    
    xv(i) = xv(i-1) + (T * fxyzRigid(xv(i-1),yv(i-1),zv(i-1),'x') );
    yv(i) = yv(i-1) + (T * fxyzRigid(xv(i-1),yv(i-1),zv(i-1),'y') );
    zv(i) = zv(i-1) + (T * fxyzRigid(xv(i-1),yv(i-1),zv(i-1),'z') );
    
end

euler = [xv; yv; zv];
euler = euler.';

end