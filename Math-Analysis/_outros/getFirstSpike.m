function first = getFirstSpike(Time_Vector)

for i=1:length(Time_Vector)

    if Time_Vector(i) ~= 0

        first = i;

        return;

    end

end


end

