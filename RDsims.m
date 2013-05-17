% Author: Andrew Harms
% Date:   03/31/2013
% Copyright (c) 2013 by Andrew Harms. This work is made available under
% the terms of the Creative Commons Attribution-ShareAlike 3.0 license

function y =...
    RDsims(num, Wvec, Rvec, Kvec, dvec, kvec, SNR, rllseq, input, ident)
% Simulations of Random Demodulator(RD) using RLL sequences.
% Reconstruction is done using Basis Pursuit and the 'yall1' package.
% 
% 'num' is the number of iterations to run for each set of parameters.
% 
% Use Wvec, Rvec, Kvec, dvec, and kvec for the desired parameter
% W, R, K, d, and k to run.  Note: 'dvec' and 'kvec' must be of equal
% length.
% 
% 'rllseq' specifies what type of random chipping sequence to use:
%   'general'
%   'repcode'
% 
% 'input' specifies whether to use a matched input signal distribution
% or a uniform input signal distribution
% 
% 'ident' is an identifier for the saved file

addpath('./YALL1-b6/');
addpath('./SpaRSA_2/');
savefile = ['./tempdata/tempsavefile' num2str(ident)];

verbose = 0; % Produce extra logs to the MATLAB command line

% number of iterations for each 4-tuple (K,R,W,d)
num_iter = num;

% RD Parameters
W_vec = Wvec;  % DFT (signal) size
R_vec = Rvec;  % Sampling rate
K_vec = Kvec;  % Sparsity level

% RLL sequence parameters: [d,k]-code
d_vec = dvec;   % min parameter
k_vec = kvec;   % max parameter

threshold = 10^(-4);  % Desired precision for reconstruction

% Multi-dim matrix for storing success probability values, or MSE for the
% noisy case, for (K,R,W,d)
prob_success_K_R_W_d =...
    zeros(length(K_vec),length(R_vec),length(W_vec),length(d_vec));

for var_W=1:length(W_vec)
    W = W_vec(var_W);
    
    % Create Unitary DFT matrix (size W x W)
    n=0:(W-1);  m=0:(W-1);
    F=exp(2*pi*1i*n'*m/W)/sqrt(W);  % Normalized DFT matrix
    if (W < R_vec(end))
        R_ind = sum((R_vec==W).*(1:length(R_vec)));
        R_use = R_vec(1:R_ind);
    else
        R_use = R_vec;
    end
    for var_R=1:length(R_use)
        R = R_use(var_R);
        r = (W/R);
        if verbose
          disp(['Running now with (W,R)=(' num2str(W) ',' num2str(R) ')']);
        end
        % Create summation-operation matrix (size RxW)
        H = create_H_matrix(R,W);
        K_max = min(length(K_vec),floor(W));
        for var_K=1:K_max
            K = K_vec(var_K);
            if verbose
              disp(['K = ' num2str(K)]);
            end
            for var_d=1:length(d_vec)
                d = d_vec(var_d);
                k = k_vec(var_d);

% ---------------------------------------------
    num_success = 0;
    for iter=1:num_iter

% -----Create Sensing Matrix-------------------
        % Create random sequence matrix D (size W x W)
        t = gen_rll_waveform(d,k,W,rllseq); % rll switching sequence
        D = diag(t); % diag matrix with switching sequence t

        % Create measurement matrix
        Phi = H*D*F;
        Phi = Phi/sqrt(r); % scale to orthonormal rows (for yall1)

% -----Create a K-sparse (column) vector-------
        leakage = 0;
        [s s_index] = create_Ksparse_vector(W,K,input,leakage);
        theta       = 2*pi*rand(1,K);
        s(s_index)  = exp(1i*theta); % random phase
        a           = ones(W,1); % unit norm entries
        s           = (a.*s);
% -----Measurements----------------------------
        y = Phi*s; % noiseless
        
        if SNR < Inf   % SNR is set to 'inf' for noiseless case
        sigPower = y'*y;
        snr      = 10^(SNR/10);
        noiseVar = sigPower./(snr*R);
        noise    = (randn(size(y)) + 1i.*randn(size(y))).*sqrt(noiseVar);
        y        = y + noise; % noisy measurements
        end

        if SNR < Inf
% -----Lasso Reconstruction--------------------
        a   = 0.9;
        tau = (1 + a)*sqrt(2*noiseVar*log(W)); % lasso parameter
        [s_hat s_hat_debias] = SpaRSA(y,Phi,tau, 'Verbose', 0,...
            'Debias', 1, 'Initialization', 2, 'StopCriterion', 1,...
            'ToleranceA', 0.01, 'MaxiterA', 100000);
        
        numNonZeroComp = sum(s_hat_debias' ~= 0);
        if verbose
            disp(['SigPower = ' num2str(sigPower) ' - EstPower = ' ...
                num2str(s_hat_debias'*s_hat_debias)...
                ' - Num NonZero Comp = ' num2str(numNonZeroComp)]);
        end
% -----Calculate error-------------------------
        if size(s_hat_debias) == size(s)
            estimate = s_hat_debias;
        else
            estimate = s_hat;
        end
        MSE    = sum((abs(estimate-s)).^2)/W;
        num_success = num_success + MSE;
        
        else % if SNR < Inf
% -----BP Reconstruction-----------------------
        opts         = struct();
        opts.tol     = 1*10^(-6);
        opts.nonorth = 0;
        s_hat        = yall1(Phi,y,opts);
% -----Calculate error-------------------------
        MSE         = max(abs(s_hat-s)); % absolute max error
        num_success = num_success + (MSE < threshold); % accuracy number
        end
    end
% ---------------------------------------------
                prob_success_K_R_W_d(var_K,var_R,var_W,var_d) =...
                    num_success/num_iter;
                % save(savefile, 'prob_success_K_R_W_d');
            end
        end
    end
end
y = prob_success_K_R_W_d;

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R TV F] = RLLspectrum(d,k,p,maxM)
% Calulate the autocorrelation (R), total variation (TV), and spectrum (F)
% of the (d,k) rll sequence with transition probability 'p'.  MaxM tells
% how much resolution is desired in the spectrum.

extraM = max(0,maxM-63);

m  = 1:min(63,maxM);
R  = zeros(1,length(m));
TV = zeros(1,length(m));

% p = [0 0.5 0.5 0.5 0.5 0.5 1];
% [P pie] = create_P(p);
[P pie] = rll_create_P(d,k,p);
[a,b]   = rll_ab_vecs(d,k,pie);
a1 = [ones(1,length(a)/2) zeros(1,length(a)/2)];
a2 = [zeros(1,length(a)/2) ones(1,length(a)/2)];
b1 = a1/21;

% Calculate the correlation R_x(m) = aP^(m)b'
for i=m
    R(i)  = a*P^(i)*b';                         % autocorrelation
    TV(i) = abs(a1*P^(i)*b1' - a2*P^(i)*b1');   % total variation
end
% Spectrum of the rll waveform
F = fft([1 R zeros(1,extraM) R(end:-1:1)])/(d+1);

% Plot spectrum if desired
plot_figures = 0;
if plot_figures
  f = (0:(length(F)-1))/length(F)-0.5;
  figure; plot(f,abs(abs((fftshift(F)))-1));
  axis([f(1) f(end) 0 2]);
  
  figure; plot(f,((abs(fftshift(F)))), 'Linewidth', 4); grid on;
  % title('Spectrum of RLL sequence');
  xlabel('Frequency Index: \omega / W', 'Fontsize', 24);
  ylabel('Amplitude', 'Fontsize', 24);
  axis([f(1) f(end) 0 2]);
  
  Frcs = 1/2 + cos(2*pi*f+pi)/2;
  figure; plot(f,((abs(fftshift(Frcs)))), 'Linewidth', 4); grid on;
  % title('Spectrum of RLL sequence');
  xlabel('Frequency Index: \omega / W', 'Fontsize', 24);
  ylabel('Amplitude', 'Fontsize', 24);
  axis([f(1) f(end) 0 2]);
end % if plot_figures

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pie = rll_create_pie(P)
% Create stationary prob distribution for RLL Markov chain from
% transition probability matrix 'P'
[V,D] = eig(P');
pie   = V(:,1)';        % stationary distribution
pie   = pie./sum(pie);  % Normalize
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
