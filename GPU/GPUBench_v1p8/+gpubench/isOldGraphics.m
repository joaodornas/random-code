function out = isOldGraphics()
%ISOLDGRAPHICS Determine whether we are using MATLABs old graphics system

%   Copyright 2014 The MathWorks, Inc.

try
    out = verLessThan('matlab','8.4.0');
catch err %#ok<NASGU>
    % If we couldn't even call the test, assume old graphics
    out = true;
end

end

