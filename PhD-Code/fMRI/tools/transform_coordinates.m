
function [x,y,z] = transform_coordinates(matrix,threshold)

    [I,J] = ind2sub(size(matrix),find(matrix < threshold));

    x = zeros(length(I),1);
    y = zeros(length(I),1);
    z = zeros(length(I),1);

    for i=1:length(I)

        x(i) = I(i);

        if J(i) <= size(matrix,2)

            y(i) = J(i);
            z(i) = 1;

        else

            rest = mod(J(i),size(matrix,2));

            if rest == 0
                
                y(i) = size(matrix,2);
                z(i) = J(i)/size(matrix,2);
                
            else
                
                quo = (J(i) - rest)/size(matrix,2);

                z(i) = quo + 1;
                y(i) = rest;
                
            end

        end

    end

end