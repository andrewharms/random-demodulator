% Author: Andrew Harms
% Date:   03/31/2013
% Copyright (c) 2013 by Andrew Harms. This work is made available under
% the terms of the Creative Commons Attribution-ShareAlike 3.0 license

function [a,b] = rll_ab_vecs(d,k,pie)
% Create RLL a,b vectors (weighted sum of symbols transmitted on
% arriving at state i or departing state j).  RLL parametes are (d,k) and
% 'pie' is the stationary distribution

if (k==inf)
    b = [ones(1,d+1) -ones(1,d+1)];
    a = b*diag(pie);
else
    b = [ones(1,k+1) -ones(1,k+1)];
    a = b*diag(pie);
end % if
end % function
