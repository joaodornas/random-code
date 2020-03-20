function [X, Y] = removeNegativos(X, Y)



    idx_X_0 = find(X<0);
    X(idx_X_0) = [];
    Y(idx_X_0) = [];
    
    idx_Y_Inf = find(Y<0);
    X(idx_Y_Inf) = [];
    Y(idx_Y_Inf) = [];
  

end

