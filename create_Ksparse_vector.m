% Author: Andrew Harms
% Date:   03/31/2013
% Copyright (c) 2013 by Andrew Harms. This work is made available under
% the terms of the Creative Commons Attribution-ShareAlike 3.0 license

function [s Tindex] = create_Ksparse_vector(W,K,input,leakage)
% Create a vector with 'K' non-zero entries of size 'W'.  'input' specifies
% matched to an rll power spectrum or uniform.  'leakage' determines
% whether to allow for off-grid frequencies
    
    switch input
        case 'matched'
            d = 1;
            k = 20;
            p = 0.6;
        case 'uniform'
            d = 0;
            k = Inf;
            p = 0.5;
        otherwise
            error('Invalid input specified');
    end
    
    % calculate the pdf 
    overSampleFactor = 4;
    [R TV pdf] = RLLspectrum(d,k,p,overSampleFactor*W);
    pdf = pdf./sum(pdf); pdf = abs(pdf); % PDF has a small imag component
    % Generate the CDF of the input distribution using the RLL pdf
    CDF = zeros(1,length(pdf));
    CDF(1) = pdf(1);
    for index=2:length(CDF)
        CDF(index) = CDF(index-1) + pdf(index);
    end
    
    % Choose K frequencies using the CDF
    p = (1:length(CDF))./length(CDF);
    u = (randsample(length(CDF),K)-1)./length(CDF);
    T = zeros(1,K);
    for i1 = 1:K
        freq = u(i1);
        index = 1;
        while freq > CDF(index)
            index = index + 1;
        end
        T(i1) = p(index);
    end
    
    s      = zeros(1,W);
    Tindex = ceil(T.*W); % index of tones
    for i1=1:length(Tindex)
        thisToneIndex = Tindex(i1);
        thisToneValue = randn(1,1) + 1i*randn(1,1); % random ampl/phase
        if leakage
            offset = rand(1,1); % random offset
        else
            offset = 0;
        end
        sTemp = zeros(1,W);
        sTemp(mod(((1:length(sTemp)) + thisToneIndex-1), W)+1) =...
            sinc((0:length(sTemp)-1) + offset)*thisToneValue; % sinc window
        % Note: can also use other windowing functions above
        s = s + sTemp;
    end
    s = s'; % output a column vector
end % function
