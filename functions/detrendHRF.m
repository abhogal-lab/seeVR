function [detrended] = detrendHRF(data, breakPT )
global opts
%setup detrend data
detrended = zeros(size(data));

xdata = [1:1:size(data,2)];
xData = xdata';
xData(breakPT(1):breakPT(2)) = [];
for ii=1:size(data,1)
input = data(ii,:)';
input(breakPT(1):breakPT(2)) = [];


%ft = fittype( 'poly1' );
%options = fitoptions( ft );

%[fitresult, ~] = fit( xData, input, ft, options );


%detrended(ii,:) = data(ii,:) - feval(fitresult,xdata)';
P = polyfit(xData,input,1)
yfit = P(1)*xdata+P(2);
detrended(ii,:) = data(ii,:) - yfit;
end
if opts.verbose
figure;
subplot(1,2,1)
for jj=1:size(data,1); data(jj,:) = rescale(data(jj,:)); end
plot(data')
subplot(1,2,2)
for jj=1:size(detrended,1); detrended(jj,:) = rescale(detrended(jj,:)); end
plot(detrended')
end
end




