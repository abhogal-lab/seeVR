function [fftexp] = fftexponential(a,s,t)

expn = s*exp(-t/a);

fftexp = fftshift(fft(expn));
end