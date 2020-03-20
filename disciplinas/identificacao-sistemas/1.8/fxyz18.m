function k = fxyz18(x,y,z,s)


% a = 1;
% b = 1;
% c = 1;
% d = 1;
% e = 1;
% f = 1;
% g = 1;
% h = 1;

a = 0.311;
b = 0.518;
c = 1.036;
d = 0.311;
e = 0.161;
f = 4.599;
g = 2.469;
h = 0.322;

switch s
    
    case 'x'
   
        k = x - x^2 + (x*y)/(x+a);
        
    case 'y'
        
        k = -b*y + (c*x*y)/(x+d) - (y*z)/(y+e);
        
    case 'z'
        
        k = f*z^2 - (g*z^2)/(y+h);
        
end


end