function testeSerial

A = rand(3,3);

B = rand(3,3);

for i=1:5
    
   disp(strcat('Divide A por B e multiplica por:',int2str(i)));
   
   p = ( A ./ B ) .* i ;
   
   p = reshape(p.',[],1);
   
   for j=1:length(p)
       
       disp(strcat('Eu ainda estou em operacao:',int2str(p(j))));
       
   end
    
end


end

