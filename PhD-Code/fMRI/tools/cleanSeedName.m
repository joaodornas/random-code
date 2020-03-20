function new_name = cleanSeedName( name )

name = name(5:end);

if strcmp('l',name(1)); name = name(3:end); end
if strcmp('r',name(1)); name = name(3:end); end

if length(name) > 2
    
    final = str2num(name(end));
    before_final = name(end-1);

    if (final >= 1)

        if (final <=5)

            if strcmp('-',before_final)

                name = name(1:end-2);

            end

        end

    end
    
end

new_name = name;

end

