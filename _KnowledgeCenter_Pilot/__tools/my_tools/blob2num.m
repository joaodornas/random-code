
function res = blob2num(x)

    res = str2double(regexp(char(x'),'[^,]+','match')');
    
end