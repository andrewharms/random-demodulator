% Author: Andrew Harms
% Date:   03/31/2013
% Copyright (c) 2013 by Andrew Harms. This work is made available under
% the terms of the Creative Commons Attribution-ShareAlike 3.0 license

function H = create_H_matrix(R,W)
% Create summation-operation matrix (size R x W).  If W/R is not an
% integer, then split entries between rows.  W and R should be integers.

if W < R
    error('Cannot use W < R');
else
    H = zeros(R,W);
    r=0; k=0;
    for i=1:R
        index1 = floor(W/R-(1-r));
        next_r = round(W-R*index1-R*(1-r))/R;
        M      = [(1-r) ones(1,index1) next_r*ones(1,next_r>0)];
        r      = next_r;
        H(i,:) = [zeros(1,k) M zeros(1,W-k-length(M))];
        if r
            k = k + (length(M)-1);
        else
            k = k + length(M);
        end
    
    end % for
end % if W < R    
end % function
