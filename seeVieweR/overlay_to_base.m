function volOnBase = overlay_to_base(baseVol, A_base, movVol , A_mov , method)
                                     
% Put movVol (T1) into the voxel grid of baseVol
%
% baseVol  : reference data   (nyB × nxB × nzB)
% A_base   : 4×4 row-vector voxel→world from base header
% movVol   : moving   data   (nyM × nxM × nzM)
% A_mov    : 4×4 row-vector voxel→world from moving header
% method   : 'nearest' | 'linear' | 'cubic'   (default 'linear')
%
% volOnBase: resampled movVol with size(baseVol)

if nargin<5, method = 'linear'; end

if isfield(A_base,'Qfactor') && isfield(A_mov,'Qfactor')
    if A_base.Qfactor ~= A_mov.Qfactor
        movVol = flip(movVol,3);                % flip slices
        Nz     = size(movVol,3);
        Fz     = eye(4);  Fz(3,3) = -1;  Fz(4,3) = Nz-1;
        A_mov  = Fz * A_mov.Transform.T;
    else
        A_mov  = A_mov.Transform.T; 
    end
end

% --- index grids of the reference volume ---------------------------------
[Xi,Yi,Zi] = ndgrid(0:size(baseVol,1)-1, ...
                    0:size(baseVol,2)-1, ...
                    0:size(baseVol,3)-1);     % 0-based like NIfTI

pts   = [Yi(:) Xi(:) Zi(:) ones(numel(Xi),1)] * A_base.Transform.T;   % world mm

% --- convert those world mm to *moving* voxel indices --------------------
idxM  = pts / A_mov;                      % row-vector algebra
xm    = idxM(:,2);      %   columns  (X)
ym    = idxM(:,1);      %   rows     (Y)
zm    = idxM(:,3);      %   slices   (Z)

% MATLAB’s interpn is 1-based.  Add 1 to go from 0-based → 1-based coords
xm = xm + 1;   ym = ym + 1;   zm = zm + 1;

% --- build scattered sampling -------------------------------------------
[NyM,NxM,NzM] = size(movVol);
volOnBase = interpn(double(movVol), ...
                    ym, xm, zm, ...         % query points
                    method, NaN);           % NaN outside FOV

volOnBase = reshape(volOnBase, size(baseVol));
volOnBase = imrotate(volOnBase,270); % needed to add this to fix orientation issues
volOnBase = flip(volOnBase,2);
end