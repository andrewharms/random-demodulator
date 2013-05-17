% Author: Andrew Harms
% Date:   03/31/2013
% Copyright (c) 2013 by Andrew Harms. This work is made available under
% the terms of the Creative Commons Attribution-ShareAlike 3.0 license

function pie = rll_create_pie(P)
% Create stationary prob distribution for RLL Markov chain from
% transition probability matrix 'P'

[V,D] = eig(P');
pie   = V(:,1)';        % stationary distribution
pie   = pie./sum(pie);  % Normalize
end % function
