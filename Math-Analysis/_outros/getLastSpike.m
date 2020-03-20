function last = getLastSpike(Time_Vector)


for i=length(Time_Vector):-1:1

    if Time_Vector(i) ~= 0

        last = i;

        return;

    end

end


end

