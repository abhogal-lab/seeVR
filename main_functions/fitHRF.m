function maps = fitHRF(mask, data, HRF_probe, HRFidx, HRFmeta, opts)
% fitHRF4 — Best-fit HRF with parameter VALUE maps (fast + block safe + saving)
%
% INPUTS
%   mask, data, HRF_probe, HRFidx, HRFmeta, opts
%   (see previous notes)
%
% OUTPUTS
%   maps.* : index maps, value maps, core maps, metadata
%
% Files are also saved to opts.hrfdir if opts.niiwrite==1
%

warning('off');
global opts;

% ---------- defaults ----------
if ~isfield(opts,'verbose'),        opts.verbose = 0; end
if ~isfield(opts,'prewhite'),       opts.prewhite = 0; end
if ~isfield(opts,'interp_factor'),  opts.interp_factor = 1; end
if ~isfield(opts,'include_linear'), opts.include_linear = 0; end
if ~isfield(opts,'niiwrite'),       opts.niiwrite = 0; end
if ~isfield(opts,'resultsdir') || isempty(opts.resultsdir), opts.resultsdir = pwd; end
if ~isfield(opts,'block_sizeV') || isempty(opts.block_sizeV), opts.block_sizeV = 50000; end
if ~isfield(opts,'mapDatatype') || isempty(opts.mapDatatype), opts.mapDatatype = 'single'; end

opts.hrfdir = fullfile(opts.resultsdir,'HRF');
if exist(opts.hrfdir,'dir') ~= 7, mkdir(opts.hrfdir); end

% ---------- shapes & ROI ----------
[x,y,z,Tdata] = size(data);
[voxel_ts, coordinates] = grabTimeseries(double(data), logical(mask)); % V×Tdata
V = size(voxel_ts,1);

% ---------- optional prewhiten ----------
if opts.prewhite
    voxel_ts = voxel_ts.';
    parfor v = 1:V
        [voxel_ts(:,v),~,~] = prewhiten(voxel_ts(:,v));
    end
    voxel_ts = voxel_ts.';
end

% ---------- align time length ----------
Tbank = size(HRF_probe,2);
if Tbank ~= Tdata
    if Tbank == Tdata*opts.interp_factor
        vt = zeros(V, Tbank);
        parfor v = 1:V, vt(v,:) = interp(voxel_ts(v,:), opts.interp_factor); end
        voxel_ts = vt; clear vt
    else
        t0 = linspace(1, Tdata, Tdata);
        t1 = linspace(1, Tdata, Tbank);
        vt = zeros(V, Tbank);
        parfor v = 1:V, vt(v,:) = interp1(t0, voxel_ts(v,:), t1, 'linear', 'extrap'); end
        voxel_ts = vt; clear vt
    end
end

% =====================================================================
% FAST path: no trend term (include_linear == 0)
% =====================================================================
if ~opts.include_linear
    X  = HRF_probe;            % K×T
    K  = size(X,1);

    % Z-score HRF bank once
    mx = mean(X,2);
    sx = std(X,0,2) + eps;
    Xz = (X - mx) .* (1./sx);  % K×T

    % Allocate outputs
    r2_best = -inf(1,V);
    Ibest   = ones(1,V);
    beta    = zeros(1,V);

    bs = opts.block_sizeV;
    for a = 1:bs:V
        b  = min(V, a+bs-1);
        Nb = b - a + 1;

        Y  = voxel_ts(a:b,:);              % Nb×T
        my = mean(Y,2);
        sy = std(Y,0,2) + eps;
        Yz = (Y - my) ./ sy;               % Nb×T

        % All correlations (K×Nb)
        corrKV_blk = (Xz * Yz.') / (size(X,2)-1);

        % Winners
        [r2_blk, i_blk] = max(corrKV_blk.^2, [], 1);        % 1×Nb

        % Slope: beta = corr * std(y) / std(x)
        corr_sel = corrKV_blk(sub2ind(size(corrKV_blk), i_blk, 1:Nb)).'; % Nb×1
        beta_blk = (corr_sel.').*( (sy.')./(sx(i_blk).') );               % 1×Nb

        % Updates
        upd = r2_blk > r2_best(a:b);
        r2_slice = r2_best(a:b); r2_slice(upd) = r2_blk(upd);  r2_best(a:b) = r2_slice;
        i_slice  = Ibest(a:b);   i_slice(upd)  = i_blk(upd);   Ibest(a:b)   = i_slice;
        b_slice  = beta(a:b);    b_slice(upd)  = beta_blk(upd); beta(a:b)    = b_slice;
    end

% =====================================================================
% Trend path: include linear term
% =====================================================================
else
    if opts.verbose, disp('fitHRF4: include_linear=1 (trend term)'); end
    K = size(HRF_probe,1);
    L = rescale(legendreP(1, 1:1:size(HRF_probe,2)), -1, 1).';  % T×1

    % Cache pseudoinverse per HRF
    Pinv = cell(K,1);
    for k = 1:K
        xk = HRF_probe(k,:).';
        Dk = [ones(size(xk)), L, xk];     % T×3
        Pinv{k} = pinv(Dk, 1e-10);        % 3×T
    end

    r2_best = -inf(1,V);
    Ibest   = ones(1,V);
    beta    = zeros(1,V);

    bs = opts.block_sizeV;
    for a = 1:bs:V
        b  = min(V, a+bs-1);
        Yt = voxel_ts(a:b,:).';           % T×Nb
        Nb = size(Yt,2);

        Ym  = mean(Yt,1);
        SST = sum( (Yt - ones(size(Yt,1),1)*Ym).^2, 1 );

        r2_blk_best = -inf(1,Nb);
        i_blk_best  = ones(1,Nb);
        beta_blk    = zeros(1,Nb);

        for k = 1:K
            B    = Pinv{k} * Yt;                   % 3×Nb
            xk   = HRF_probe(k,:).';
            Dk   = [ones(size(xk)), L, xk];        % T×3
            Yhat = Dk * B;                         % T×Nb

            SSE  = sum( (Yt - Yhat).^2, 1 );
            R2k  = 1 - (SSE ./ (SST + eps));

            upd = R2k > r2_blk_best;
            r2_blk_best(upd) = R2k(upd);
            i_blk_best(upd)  = k;
            beta_blk(upd)    = B(3,upd);
        end

        r2_best(a:b) = r2_blk_best;
        Ibest(a:b)   = i_blk_best;
        beta(a:b)    = beta_blk;
    end
end

% ---------- assemble core volumes ----------
HRF_index = zeros(1,x*y*z); HRF_index(1,coordinates) = Ibest; HRF_index = reshape(HRF_index,[x y z]);
r2_vol    = zeros(1,x*y*z); r2_vol(1,coordinates)    = r2_best; r2_vol  = reshape(r2_vol,[x y z]);
beta_vol  = zeros(1,x*y*z); beta_vol(1,coordinates)  = beta;    beta_vol= reshape(beta_vol,[x y z]);

% ---------- parameter INDEX & VALUE maps ----------
best_idx = HRFidx(Ibest,:);     % V×(2 or 3)

onset_idx = best_idx(:,1);
disp_idx  = best_idx(:,2);
if size(best_idx,2) >= 3
    under_idx = best_idx(:,3);
else
    under_idx = ones(size(Ibest)); % dummy
end

% Values
if isfield(HRFmeta,'use_onset') && HRFmeta.use_onset
    onset_val = zeros(size(onset_idx));
    nz = onset_idx > 0;
    if any(nz)
        onset_vals_grid = HRFmeta.onset_values(:);
        onset_val(nz)   = onset_vals_grid(onset_idx(nz));
    end
else
    onset_val = zeros(size(onset_idx));
end

disp_vals_grid = HRFmeta.disp_values(:);
disp_val       = disp_vals_grid(disp_idx);

if strcmpi(HRFmeta.model,'gamma') && isfield(HRFmeta,'under_values')
    under_vals_grid = HRFmeta.under_values(:);
    under_val       = under_vals_grid(under_idx);
else
    under_val = zeros(size(under_idx));
end

% Scatter to volumes
onset_idx_vol = zeros(1,x*y*z); onset_idx_vol(1,coordinates) = onset_idx; onset_idx_vol = reshape(onset_idx_vol,[x y z]);
disp_idx_vol  = zeros(1,x*y*z); disp_idx_vol(1,coordinates)  = disp_idx;  disp_idx_vol  = reshape(disp_idx_vol,[x y z]);
under_idx_vol = zeros(1,x*y*z); under_idx_vol(1,coordinates) = under_idx; under_idx_vol = reshape(under_idx_vol,[x y z]);

onset_val_vol = zeros(1,x*y*z); onset_val_vol(1,coordinates) = onset_val; onset_val_vol = reshape(onset_val_vol,[x y z]);
disp_val_vol  = zeros(1,x*y*z); disp_val_vol(1,coordinates)  = disp_val;  disp_val_vol  = reshape(disp_val_vol,[x y z]);
under_val_vol = zeros(1,x*y*z); under_val_vol(1,coordinates) = under_val; under_val_vol = reshape(under_val_vol,[x y z]);

% ---------- save out ----------
isExp = isfield(HRFmeta,'model') && strcmpi(HRFmeta.model,'exp');
if opts.niiwrite
    cd(opts.hrfdir);
    if isExp
        niftiwrite(cast(HRF_index,opts.mapDatatype),'EXP_HRF_index',opts.info.map);
        niftiwrite(cast(r2_vol,   opts.mapDatatype),'EXP_HRF_r2',   opts.info.map);
        niftiwrite(cast(beta_vol, opts.mapDatatype),'EXP_HRF_beta', opts.info.map);
        niftiwrite(cast(disp_val_vol,  opts.mapDatatype),'EXP_tau_value',    opts.info.map);
        if isfield(HRFmeta,'use_onset') && HRFmeta.use_onset
            niftiwrite(cast(onset_val_vol, opts.mapDatatype),'EXP_delay_value', opts.info.map);
        end
    else
        niftiwrite(cast(HRF_index,opts.mapDatatype),'GAMMA_HRF_index',opts.info.map);
        niftiwrite(cast(r2_vol,   opts.mapDatatype),'GAMMA_HRF_r2',   opts.info.map);
        niftiwrite(cast(beta_vol, opts.mapDatatype),'GAMMA_HRF_beta', opts.info.map);

        niftiwrite(cast(onset_idx_vol,opts.mapDatatype),'GAMMA_onset_idx', opts.info.map);
        niftiwrite(cast(disp_idx_vol, opts.mapDatatype),'GAMMA_disp_idx',  opts.info.map);
        niftiwrite(cast(under_idx_vol,opts.mapDatatype),'GAMMA_under_idx', opts.info.map);

        niftiwrite(cast(onset_val_vol,opts.mapDatatype),'GAMMA_onset_value', opts.info.map);
        niftiwrite(cast(disp_val_vol, opts.mapDatatype),'GAMMA_disp_value',  opts.info.map);
        niftiwrite(cast(under_val_vol,opts.mapDatatype),'GAMMA_under_value', opts.info.map);
    end
else
    if isExp
        saveImageData(HRF_index,opts.headers.map,opts.hrfdir,'EXP_HRF_index.nii.gz',64);
        saveImageData(r2_vol,   opts.headers.map,opts.hrfdir,'EXP_HRF_r2.nii.gz',64);
        saveImageData(beta_vol, opts.headers.map,opts.hrfdir,'EXP_HRF_beta.nii.gz',64);
        saveImageData(disp_val_vol, opts.headers.map,opts.hrfdir,'EXP_tau_value.nii.gz',64);
        if isfield(HRFmeta,'use_onset') && HRFmeta.use_onset
            saveImageData(onset_val_vol,opts.headers.map,opts.hrfdir,'EXP_delay_value.nii.gz',64);
        end
    else
        saveImageData(HRF_index,opts.headers.map,opts.hrfdir,'GAMMA_HRF_index.nii.gz',64);
        saveImageData(r2_vol,   opts.headers.map,opts.hrfdir,'GAMMA_HRF_r2.nii.gz',64);
        saveImageData(beta_vol, opts.headers.map,opts.hrfdir,'GAMMA_HRF_beta.nii.gz',64);
        saveImageData(onset_idx_vol,opts.headers.map,opts.hrfdir,'GAMMA_onset_idx.nii.gz',64);
        saveImageData(disp_idx_vol, opts.headers.map,opts.hrfdir,'GAMMA_disp_idx.nii.gz',64);
        saveImageData(under_idx_vol,opts.headers.map,opts.hrfdir,'GAMMA_under_idx.nii.gz',64);
        saveImageData(onset_val_vol,opts.headers.map,opts.hrfdir,'GAMMA_onset_value.nii.gz',64);
        saveImageData(disp_val_vol, opts.headers.map,opts.hrfdir,'GAMMA_disp_value.nii.gz',64);
        saveImageData(under_val_vol,opts.headers.map,opts.hrfdir,'GAMMA_under_value.nii.gz',64);
    end
end

% ---------- pack outputs ----------
maps = struct();
maps.HRF_index = HRF_index;
maps.r2        = r2_vol;
maps.beta      = beta_vol;
maps.onset_idx = onset_idx_vol;
maps.disp_idx  = disp_idx_vol;
maps.under_idx = under_idx_vol;
maps.onset_val = onset_val_vol;
maps.disp_val  = disp_val_vol;
maps.under_val = under_val_vol;
maps.meta      = HRFmeta;
end
