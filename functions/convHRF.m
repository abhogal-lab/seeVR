%Copyright Alex A. Bhogal, 7/15/2021, University Medical Center Utrecht, 
%a.bhogal@umcutrecht.nl
%The seeVR toolbox is software, licensed under the Creative Commons 
%Attribution-NonCommercial-ShareAlike 4.0 International Public License
%By using seeVR and associated scripts you agree to the license conditions
%that can be reviewed at:
%https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode
%These tools are for research purposes and are not intended for
%commercial purposes. 

function [HRF,HRF_probe] = convHRF(HRFprobe, opts)

if isfield(opts,'rratio'); else; opts.rratio = 1000; end %seems to affect the vertical spread of HRF (use large value to limit)

input_probe = HRFprobe;
pad = 2^nextpow2(length(HRFprobe));
HRFprobe = zeros(length(input_probe)+2*pad,1);
HRFprobe(1:pad) = input_probe(1,1);
HRFprobe(pad+1:end-pad,1) = input_probe;
HRFprobe(end-pad+1:end,1) = input_probe(1,end)
mm=0
xdata = [1:1:length(HRFprobe)];
clear HRF


for ii = opts.onset %onset
    for jj =  opts.disp %dispersion
        for kk = opts.under %undershoot
            mm = mm+1;
            V = [ii jj kk];          
            t = [0:1:length(xdata)-1];
            h_final = zeros(length(t),1);
            alpha_1 = V(1);
            beta_1 = V(2);
            alpha_2 = alpha_1;
            beta_2 = V(3);
            
            h = ((t.^(alpha_1 - 1)).*(((beta_1).^(-alpha_1)).*(exp((1./(-beta_1)).*t))))./(gamma(alpha_1));
            h2 = ((t.^(alpha_2 - 1)).*(((beta_2).^(-alpha_2)).*(exp((1./(-beta_2)).*t))))./(opts.rratio.*gamma(alpha_2));
            
            h_final = h - h2;
            h_final = h_final./trapz(h_final);
            
            HRF(mm,:) = h_final;
            
        end
    end
end

% convolve HRF with probe
HRF_probe = zeros([mm+1 length(HRFprobe)]);
HRF_probe(1,:) = HRFprobe/trapz(HRFprobe);
newprobe_freq = fftshift(fft(HRFprobe)/trapz(HRFprobe));

for ii=1:size(HRF,1)
    HRF_freq(ii,:) = fftshift(fft(HRF(ii,:)/trapz(HRF(ii,:))));
    tmp = real(ifft(ifftshift(HRF_freq(ii,end:-1:1).*newprobe_freq')));
    HRF_probe(ii+1,:) = (tmp/trapz(tmp)); %normalize HRF
    HRF_probe(ii+1,:) =   HRF_probe(ii+1,end:-1:1);
end

HRF_probe(:,1:pad) = []; HRF_probe(:,end-pad+1:end) = [];  %remove padding
for ii=1:size(HRF_probe,1); HRF_probe(ii,:) = rescale(HRF_probe(ii,:)); end
if opts.plot
figure; 
customMap = viridis(size(HRF,1)); customMap = flip(customMap,1);
for ii=1:size(HRF,1)
subplot(1,2,1); plot(t,HRF(ii,:)', 'Color', customMap(ii,:)); ylim([-0.2 1]); xlim([0 50]); title('HRF'); hold on;
subplot(1,2,2); plot(HRF_probe(ii+1,:)', 'Color', customMap(ii,:)); hold on;
end
plot( HRF_probe(1,:), 'k', 'LineWidth', 2);  title('HRF probe');
if isfield(opts,'figdir')
saveas(gcf,[opts.figdir,'HRF_fnc.fig']);
else
saveas(gcf,[pwd,'\','HRF_fnc.fig']);  
end
end