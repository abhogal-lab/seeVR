function [newprobe] = optimizeRegressor(probe, data, refmask, opts)

if ~isfield(opts,'rescale_probe'); opts.rescale_probe = 1; end       %rescaling may be helful for initial refinement run
if ~isfield(opts,'corr_thresh'); opts.corr_thresh = 0.7; end         %threshold by which to accept correlated voxels during the refinement of the BOLD regressor
if ~isfield(opts,'lowerthresh'); opts.lowerthresh = -2; end
if ~isfield(opts,'upperthresh'); opts.upperthresh = 2; end
if ~isfield(opts,'pca_analysis'); opts.pca_analysis = 1; end
if ~isfield(opts,'RMSEthresh'); opts.RMSEthresh = 0.05; end
if ~isfield(opts,'explainedVariance'); opts.explainedVariance = 0.75; end

probe = probe(:);

[voxel_ts,~] = grabTimeseries(data, refmask);

disp('Creating optimized regressor')

keep_probes = [];
%keeps looping until RSME between probe is very small = therefore, reach
%convergence!
roundprobe=0;
stop=0;
while stop==0
    tic
    roundprobe=roundprobe+1;

    if roundprobe == 1 %rescaling probe may be helpful for initial iteration
        if opts.rescale_probe
            newprobe = rescale(probe);
        else
            newprobe = probe;
        end
    end

    disp(['Correlating TS voxelwise to probe ',num2str(roundprobe)]);disp(' . ');disp(' . ');disp(' . ');
    % voxelwise corrrelation with probe

    [~,lags] = xcorr(newprobe,voxel_ts(1,:)','coeff');
    %matrix cross correlation with probe

    a=voxel_ts;
    a2=mat2cell(a,ones(1,size(a,1)),size(a,2)); %turn into cells so treat all rows independently

    % the order here matters because of the SIGN of the lag
    % this way we slide the probe to the ts
    b2=cellfun(@(a2) xcorr(a2,newprobe,'none'),a2,'uniformoutput',false);
    corr_probe=cell2mat(b2); %reformat into correlation matrix
    corr_probe = corr_probe' ;

    %remove low and high lags ignoring unreasonable and long lag times
    idx = find(lags<=opts.lowerthresh | lags>=opts.upperthresh);
    lags(idx)=[];
    corr_probe(idx,:)=[];      clear idx

    % estimate peak
    % delay and correlation over specified time lag range with fit parameters
    maxcorr = max(corr_probe,[],1);
    tmp = repmat(maxcorr,length(corr_probe(:,1)),1);
    yes = zeros(size(corr_probe));
    yes(find(tmp==corr_probe))=1; clear tmp
    yes(yes==0)=NaN;
    tmp = repmat(lags',1,length(corr_probe(1,:)));
    ts_lag = nanmean(tmp.*yes); %clean things up

    checkcorr=maxcorr;
    %creating a filter to remove the timeseries before the PCA
    checkcorr(maxcorr>=opts.corr_thresh)=0;  % if good, put 0 because will be kept
    checkcorr(maxcorr<opts.corr_thresh)=1;  % to be removed

    checklag=ts_lag;
    checklag(ts_lag>=opts.lowerthresh & ts_lag<=opts.upperthresh)=0;
    checklag(ts_lag<opts.lowerthresh | ts_lag>opts.upperthresh)=1;

    resultant_filter = zeros(size(checklag));
    resultant_filter(checklag==1 | checkcorr==1)=1;
    remove = find(resultant_filter==1);

    % removing timeseries and lags outside of range
    voxel_ts_filt=voxel_ts';
    voxel_ts_filt(:,remove)=[];
    ts_lag_filt=ts_lag;
    ts_lag_filt(remove)=[];

    % re-aligning all timeseries based on lag for PCA

    new_shifted_ts = zeros(length(ts_lag_filt), length(probe));
    for hh = 1:length(ts_lag_filt)
        new_shifted_ts(hh,:) = circshift(voxel_ts_filt(:,hh),ts_lag_filt(hh));
        if ts_lag_filt(hh) > 0
            new_shifted_ts(hh,1:ts_lag_filt(hh)) = voxel_ts_filt(1,1);
        elseif ts_lag_filt(hh) < 0
            new_shifted_ts(hh,end-ts_lag_filt(hh):end) = voxel_ts_filt(1,end);
        else
            new_shifted_ts(hh,:) = voxel_ts_filt(:,hh);
        end
    end

    disp('finished correlation, estimating new probe')

    if opts.pca_analysis
        try
            %%% PCA on significant shifted timeseries to get single timecourse that explains
            %the highest shared variance
            % score = how much pattern from each component is in the voxel
            % coeff = actual component pattern

            [coeff,score,latent,~,explained] = pca(new_shifted_ts,'VariableWeights','variance','Row','Complete'); % complete option ignores NaN

            %loop until have enough components to explain 85% of the
            %variability
            total_components = 3;jj=total_components-1;  % minimum 3 components
            while jj < length(explained)
                jj=jj+1;
                checkk = sum(explained(1:jj));
                if checkk >= opts.explainedVariance % if first 3 component covers more than 90, continue
                    continue
                else
                    total_components = total_components+1; % see explained, usually explains ~ 75-85% with first 10 components
                end
            end

            try
                pca_components = score(:,1:total_components)*coeff(:,1:total_components)';
            catch
                disp('The correlation parameter for the optimized regressor is set to high')
                disp('i.e. there are not sufficient highly correlated voxels: set opts.corrthresh to a lower value')
                disp(['the opts.corr_thresh parameter is currently set at: ',num2str(opts.corr_thresh)])
                error('exiting function')
                return
            end
            pca_components = pca_components(~any(isnan(pca_components),2),:); %removing rows with NaN

            % components can be flipped so doing a pearson correlation with the
            % probe gas trace; if flipped, will have a - correlation so
            % just multiply by -1
            outputpca = nan(size(pca_components));
            clear pearsR
            pearsR = zeros([1 length(pca_components(:,1))]);
            parfor dd = 1:length(pca_components(:,1))
                pearsR(dd) = corr(newprobe,pca_components(dd,:)');

                if pearsR(dd) > 0
                    outputpca(dd,:) = 1 * pca_components(dd,:);
                else
                    outputpca(dd,:) = -1 * pca_components(dd,:);
                end

            end

            % weigthed average of refined regressor
            newprobe=[];
            newprobe = nanmean(outputpca,1)';
            toc
        catch
            disp('probably cant do PCA due to missing toolbox')
        end
    else
        newprobe=[];
        newprobe = (nanmean(new_shifted_ts,1))';
    end

    % averaging all targeted components for new probe
    keep_probes(roundprobe,:) = newprobe;
    clear max_corr_voxel   pca_components  outputpca   new_shifted_ts pearsR

    rmse=[];
    if roundprobe==1
        rmse = immse(probe(:),newprobe(:))
    else
        rmse = immse(newprobe',(keep_probes(roundprobe-1,:)))
    end

    if rmse>0 & rmse<opts.RMSEthresh
        stop = 1;
    elseif roundprobe==10 %maximum number of iterations to find probe
        stop = 1;
    end

end
figure(15);
for ii=1:size(keep_probes,1); legendInfo{ii} = ['probe ',int2str(ii)]; end
subplot(3,1,1);
plot(probe,'LineWidth',2); title('input probe' )
subplot(3,1,2);
plot(keep_probes','LineWidth',2); title('probe iterations');
legend(legendInfo)
subplot(3,1,3);
plot(newprobe,'LineWidth',2); title('optimized probe');
saveas(gcf,fullfile(opts.resultsdir,'all_probes.fig'));
%save the final probe
save(fullfile(opts.resultsdir,'final_probe.mat'), 'newprobe');

disp('Finished creating optimized regressor')

end

