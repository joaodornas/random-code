function [width] = fullWidthAtHalfHeight(x_vector,y_vector,scaleMode)

%[width] = fullWidthAtHalfHeight(vector,scaleMode)
%calculates full width at half height of a function whose y axis is given
%in y_vector and x axis in x_vector. scaleMode is a string input that
%determines if scale should be linear ('linear'), octave (log2 - 'octave')
%or log ('log').
%
%by Lucas Pinto in May 2009

y_vector=fixdec(y_vector,1);
maxValue = max(y_vector);
i_maxValue = find(y_vector==maxValue,1,'first');
i_half1 = find(y_vector==fixdec(maxValue/2,1),1,'first');
i_half2 = find(y_vector==fixdec(maxValue/2,1),1,'last');

if isempty(i_half1)==1 || isempty(i_half2)==1
    width=nan;
elseif i_half1>i_maxValue || i_half2<i_maxValue
    width=nan;
else

    switch scaleMode
        case 'linear'
            width=x_vector(i_half2)-x_vector(i_half1);
        case 'octave'
            width=log2(x_vector(i_half2))-log2(x_vector(i_half1));
        case 'log'
            width=log(x_vector(i_half2))-log(x_vector(i_half1));
    end

end