% Author: Andrew Harms
% Date:   03/31/2013
% Copyright (c) 2013 by Andrew Harms. This work is made available under
% the terms of the Creative Commons Attribution-ShareAlike 3.0 license

function y = gen_rll_seq(d,k,N)
% y = gen_rll_seq(d,k,N)
% 
% Generate an RLL (d,k) sequence on {0,1} of length 'N'
% The output is the coded sequence with random data as input

if (k < d)
    error('Invalid arguments: require d<=k');
end
y = zeros(1,N);         % Output sequence
x = round(rand(1,N));   % Random input sequence

s = 1;  % index
t = 1;  % stopping criterion
m = 0;  % track consecutive zeros
while t
    if (x(s)) % start with d zeros
        y(s) = 1;
        y((s+1):(s+d)) = zeros(1,d);
        s = s+d+1;
        m = d;
    elseif (~x(s) && (m<k) ) % generate a zero with prob. 1/2
        y(s) = 0;
        s = s+1;
        m = m+1;
    else                    % generate a one with prob. 1/2, add d zeros
        y(s) = 1;
        y((s+1):(s+d)) = zeros(1,d);
        s = s+d+1;
        m = d;
    end
    
    if (s > N) % stopping criterion at length N
        t = 0;
        y = y(1:N);
    end
end % while
end % function
