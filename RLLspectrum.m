% Author: Andrew Harms
% Date:   03/31/2013
% Copyright (c) 2013 by Andrew Harms. This work is made available under
% the terms of the Creative Commons Attribution-ShareAlike 3.0 license

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
