function [fftexp] = fftexponential(a,t)

expn = exp(-t/a);

fftexp = fftshift(fft(expn));
end