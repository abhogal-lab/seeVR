
Fs = 1/opts.TR;              
Ts = opts.TR;
fc = opts.dyn;    
V = 0.1;                 % Window Stop Time
ActTD = 0.008;           % Actual Time Delay
t = opts.xdata;
x = probe;
y =  orig_voxel_ts;
%%Spectra of Input Signals
%Correlation Length
corrLength=2*length(x)-1;
%Input Signal-1 Spectra
[~,lags] = xcorr(x,y(1,:));
X=fft(x, corrLength);
%Input Signal-2 Spectra
%Y = zeros([corrLength length(y)]);

Y=fft(y', corrLength);
%%Cross Correlation in Frequency Domain
%Hadamard Product
%the cross-correlation between two signals is equal to the product of 
%fourier transform of one signal multiplied by complex conjugate of fourier
%transform of another signal. After doing this, when we take the ifft of 
%the product signal, we get a peak which indicates the shift between two signals
Z=X.*(conj(Y));
%After doing this, when we take the ifft of the product signal, we get a 
%peak which indicates the shift between two signals
z=real(ifftshift(ifft(Z)));
for ii=1:20; figure(11); hold on; plot(rescale(z(:,1000+ii))); end
%Time Axis
[~,idx] = max(abs(z),[],1);

TD=lags(idx);     % Calculated Time Delay

[xx yy zz] = size(mask);
lgm = zeros([1 xx*yy*zz]);
lgm(1,coordinates) = TD;
lag_map = reshape(lgm,size(mask));
figure; imagesc(lag_map(:,:,20))

%%%


