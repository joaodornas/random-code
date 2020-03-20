function [X, Y] = removeInfZeros(X, Y)



    idx_X_0 = find(X==0);
    X(idx_X_0) = [];
    Y(idx_X_0) = [];
    
    idx_X_Inf = find(X==Inf);
    X(idx_X_Inf) = [];
    Y(idx_X_Inf) = [];
    
    idx_X_Minus_Inf = find(X==-Inf);
    X(idx_X_Minus_Inf) = [];
    Y(idx_X_Minus_Inf) = [];
    
    idx_Y_0 = find(Y==0);
    X(idx_Y_0) = [];
    Y(idx_Y_0) = [];
    
    idx_Y_Inf = find(Y==Inf);
    X(idx_Y_Inf) = [];
    Y(idx_Y_Inf) = [];
    
    idx_Y_Minus_Inf = find(Y==-Inf);
    X(idx_Y_Minus_Inf) = [];
    Y(idx_Y_Minus_Inf) = [];


end

