% Author: Andrew Harms
% Date:   03/31/2013
% Copyright (c) 2013 by Andrew Harms. This work is made available under
% the terms of the Creative Commons Attribution-ShareAlike 3.0 license

function [P pie] = rll_create_P(d,k,p)
% Create RLL transition matrix and stationary distribution

if (k==inf)
    P1 = eye(2*d+1);
    P = [[zeros(2*d+1,1) P1]; zeros(1,2*d+2)];
    P(d+1,d+1)     = p; P(d+1,d+2) = 1-p;
    P(2*d+2,2*d+2) = p; P(2*d+2,1) = 1-p;
    
    x = (1-p)./(2+2*d*(1-p));
    
    pie = x.*ones(1,2*d+2);
    pie(d+1) = pie(d+1)./(1-p);
    pie(end) = pie(end)./(1-p);
    
else
    p_trans = p.*ones(1,k-d);
    
    p0 = [ones(1,d) p_trans];   % Transition probabilities for symbol zero
    p1 = ones(1,k) - p0;        % Transition probabilities for symbol one
    
    % Transition matrix p_ij = P{i->j|i}
    P11 = [ [zeros(1,length(p1))' diag(p0)]; [0 zeros(1,k)] ];
    P01 = zeros(k+1); P01(k+1,1) = 1; P01(1:k,1) = p1';
    
    P = [P11 P01; P01 P11];
    
    pie = rll_create_pie(P);
end
end % function
