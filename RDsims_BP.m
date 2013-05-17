% Author: Andrew Harms
% Date:   03/31/2013
% Copyright (c) 2013 by Andrew Harms. This work is made available under
% the terms of the Creative Commons Attribution-ShareAlike 3.0 license

function y =...
    RDsims_BP(num, Wvec, Rvec, Kvec, dvec, kvec, rllseq, input, ident)
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

addpath('./functions/YALL-b6/');
savefile = ['./tempdata/tempsavefile' num2str(ident)];

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

% Multi-dim matrix for storing success probability values for (K,R,W,d)
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
        disp(['Running now with (W,R)=(' num2str(W) ',' num2str(R) ')']);
        % Create summation-operation matrix (size RxW)
        H = create_H_matrix(R,W);
        K_max = min(length(K_vec),floor(W));
        for var_K=1:K_max
            K = K_vec(var_K);
            disp(['K = ' num2str(K)]);
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

% -----BP Reconstruction-----------------------
        opts         = struct();
        opts.tol     = 1*10^(-6);
        opts.nonorth = 0;
        s_hat        = yall1(Phi,y,opts);
% -----Calculate error-------------------------
        MSE         = max(abs(s_hat-s)); % absolute max error
        num_success = num_success + (MSE < threshold); % accuracy number
    end
% ---------------------------------------------
                prob_success_K_R_W_d(var_K,var_R,var_W,var_d) =...
                    num_success/num_iter;
                save(savefile, 'prob_success_K_R_W_d');
            end
        end
    end
end
y = prob_success_K_R_W_d;
