%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Pink Noise Generation with MATLAB Implementation   %
%                                                      %
% Author: M.Sc. Eng. Hristo Zhivomirov        07/30/13 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = pinknoise_my(N,f_s,power)

% function: y = pinknoise(N) 
% N - number of samples to be returned in a row vector
% basicfre - the basic frequency of the series
% y - a row vector of pink (flicker) noise samples
% fs - the sampling frequency of the data
% power - the power at 1 Hz of the pink noise for signal has length of 1s

% The function generates a sequence of pink (flicker) noise samples. 
% In terms of power at a constant bandwidth, pink noise falls off at 3 dB/oct, i.e. 10 dB/dec. 

% difine the length of the vector and
% ensure that the M is even, this will simplify the processing
if rem(N, 2)
    M = N+1;
else
    M = N;
end

% generate white noise
x = randn(1, M);

% FFT
X = fft(x);
p2 = X/(M/f_s);
p1 = p2(1:M/2+1);
% prepare a vector with frequency indexes 
NumUniquePts = M/2 + 1;     % number of the unique fft points
n = 1:NumUniquePts;         % vector with frequency indexes 
df = f_s/M;                 % frequency increment;
f = f_s*(0:(M/2))/M;
% manipulate the left half of the spectrum so the PSD
% is proportional to the frequency by a factor of 1/f, 
% i.e. the amplitudes are proportional to 1/sqrt(f)
p1(1) = 0;  
p1(2:end) = p1(2:end)./sqrt(f(2:end));
% scale the spectrum to ensure the power at 1Hz
[~,minind] = min(abs(f-1));

ptemp = [p1 conj(p1(end-1:-1:2))];
ytemp = real(ifft(ptemp));

sca_pow = abs(p1(minind))/sqrt(power*(M/f_s)/df);
p1 = p1./sca_pow;


% prepare the right half of the spectrum - a conjugate copy of the left one,
% except the DC component and the Nyquist component - they are unique
% and reconstruct the whole spectrum
p1 = [p1 conj(p1(end-1:-1:2))];

% IFFT
y = real(ifft(p1));

ptest = abs(fft(y));

% ensure that the length of y is N
y = y(1:N);

% ensure zero mean value
y = y - mean(y);

end