function [input_freq] = fftinput(input)

%input_freq = fftshift(fft(input/trapz(input)));
input_freq = fftshift(fft(input));

end