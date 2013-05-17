% Author: Andrew Harms
% Date:   03/31/2013
% Copyright (c) 2013 by Andrew Harms. This work is made available under
% the terms of the Creative Commons Attribution-ShareAlike 3.0 license

function t = gen_rll_waveform(d,k,n,rllseq)
% Generate a modulated RLL (d,k) waveform of {+1,-1} of length 'n'
% y is the coded sequence with random data as input.
% Note: 'rllseq' specifies a "Markov generated" rll sequence or a
% repetition-coded rll sequence

switch rllseq
    case 'general'
              
    t_01 = gen_rll_seq(d,k,n);  % t_01 = [d,k] sequence of 0,1
    t    = zeros(1,n);          % t = switching sequence of +1,-1
    t(1) = 2*t_01(1)-1;         % random starting symbol
    for i2=2:length(t_01)       % diffentially encode t_01 for t
        if t_01(i2) % ones indicate transitions
            t(i2) = -sign(t(i2-1));
        else
            t(i2) = sign(t(i2-1));
        end
    end % for i2
    
    case 'repcode'
% Generate a repetition-coded waveform of {+1,-1} of length 'n' with 
% repetition rate 'd' and maximum # of repeated entries 'k'
% y is the coded sequence with random data as input
    tmp  = gen_rll_seq(0,k,n);  % [0,k] sequence of 0,1
    t_01 = zeros(1,(d+1)*n);
    t_01(1:(d+1):end) = tmp;    % ones are separated by at least 'd' zeros
    t_01 = t_01(1:n);
    t    = zeros(1,n);          % t = switching sequence of +1,-1
    t(1) = 2*t_01(1)-1;         % random initial symbol
    for i2=2:length(t_01)       % diffentially encode t_01 for t
        if t_01(i2) % ones indicate transitions
            t(i2) = -t(i2-1);
        else
            t(i2) = t(i2-1);
        end
    end % for i2
    
    otherwise
        % Check for valid sequence type
        error('Invalid rllseq specified.  Use general or repcode');
end % switch rllseq
