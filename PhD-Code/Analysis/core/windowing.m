function [o] = windowing( window, n_windows, vector )

o.window = window;
o.n_windows = n_windows;
o.middle = floor(length(vector)/2) + 1;
o.side = o.middle - 1 - floor(o.window/2);
o.brick = floor(o.side / ((o.n_windows - 1)/2));
o.left_vector = o.middle:-o.brick:(o.middle - ((o.n_windows - 1)/2)*o.brick + ((n_windows - 1)/2));
o.right_vector = o.middle:o.brick:(o.middle + ((o.n_windows - 1)/2)*o.brick - ((n_windows - 1)/2));

remove_right = o.right_vector(2:end);

o.steps = [o.left_vector, remove_right];

end

