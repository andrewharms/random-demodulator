% Author: Andrew Harms
% Date:   03/31/2013
% Copyright (c) 2013 by Andrew Harms. This work is made available under
% the terms of the Creative Commons Attribution-ShareAlike 3.0 license

function y =...
    RDsims_Lasso(num, Wvec, Rvec, Kvec, dvec, kvec, SNR, rllseq, input, ident)
% Simulations of Random Demodulator(RD) using RLL sequences.
% Reconstruction is done using the SpaRSA package (essentially Lasso).
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

addpath('./SpaRSA_2/');

verbose = 0; % Produce extra logs to the MATLAB command line
savefile = ['./tempdata/tempsavefile' num2str(ident)]; % File to save data

% number of iterations for each 4-tuple (K,R,W,d)
num_iter = num;

% RD Parameters
W_vec = Wvec; % DFT (signal) size
R_vec = Rvec; % Sampling rate
K_vec = Kvec; % Sparsity level

% RLL sequence parameters: [d,k]-code
d_vec = dvec; % min parameter
k_vec = kvec; % max parameter

s_save = zeros(num_iter, max(W_vec)); % save input signals

% Multi-dim matrix for storing success probability values for (K,R,W,d)
avgMSE_K_R_W_d =...
    zeros(length(K_vec),length(R_vec),length(W_vec),length(d_vec));

for var_W=1:length(W_vec)
    W = W_vec(var_W);
    
    % Create Unitary DFT matrix (size W x W)
    n=0:(W-1);  m=0:(W-1);
    F=exp(2*pi*1i*n'*m/W)/sqrt(W);  % Normalized DFT matrix
    % Check that R < W
    if (W < R_vec(end))
        R_ind = sum((R_vec==W).*(1:length(R_vec)));
        R_use = R_vec(1:R_ind);
    else
        R_use = R_vec;
    end
    for var_R=1:length(R_use)
        R = R_use(var_R);
        r = (W/R);
        disp(['Running now with (W,R)=(' num2str(W) ',' num2str(R) ')']);
        % Create summation-operation matrix (size RxW)
        H = create_H_matrix(R,W);
        K_max = min(length(K_vec),floor(R/4));
        for var_K=1:K_max
            K = floor(K_vec(var_K));
            disp(['K = ' num2str(K)]);
            for var_d=1:length(d_vec)
                d = d_vec(var_d);
                k = k_vec(var_d);

% ---------------------------------------------
    avgMSE = 0;
    for iter=1:num_iter

% -----Create Sensing Matrix-------------------
        % Create random sequence matrix D (size W x W)
        t = gen_rll_waveform(d,k,W,rllseq); % rll waveform with +1 and -1
        D = diag(t); % diag matrix with switching sequence t

        % Create measurement matrix
        Phi = H*D*F;
        Phi = Phi/sqrt(r);      % scale to orthonormal rows

% -----Create a K-sparse (column) vector-------
        leakage = 0;
        [s s_index] = create_Ksparse_vector(W,K,input,leakage);
        s_save(iter,1:length(s)) = s;
        
% -----Measurements----------------------------
        y        = Phi*s; % noiseless measurements
        sigPower = y'*y;
        snr      = 10^(SNR/10);
        noiseVar = sigPower./(snr*R);
        noise    = (randn(size(y)) + 1i.*randn(size(y))).*sqrt(noiseVar);
        y        = y + noise; % noisy measurements

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
        avgMSE = avgMSE + MSE;
    end %for iter=...
% ---------------------------------------------
                avgMSE_K_R_W_d(var_K,var_R,var_W,var_d) =...
                    avgMSE/num_iter;
                save(savefile, 'avgMSE_K_R_W_d');
            end
        end
    end
end

y = avgMSE_K_R_W_d;
