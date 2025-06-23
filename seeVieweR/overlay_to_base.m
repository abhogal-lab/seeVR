% Copyright (C) Alex A. Bhogal, 2025, University Medical Center Utrecht
% a.bhogal@umcutrecht.nl
% <overlay_to_base: resamples image for overlay using header information
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
function volOnBase = overlay_to_base(baseVol, baseInfo, movVol, movInfo, method)
% OVERLAY_TO_BASE Resamples movVol into baseVol's voxel grid using affine transforms.
%
% baseVol : reference volume [Y X Z]
% baseInfo: niftiinfo struct for base volume
% movVol  : moving volume [Y X Z]
% movInfo : niftiinfo struct for moving volume
% method  : interpolation method (default: 'linear')

    if nargin < 5, method = 'cubic'; end

    % Extract affine matrices
    A_base = getAffineMatrix(baseInfo);
    A_mov  = getAffineMatrix(movInfo);

    % Generate voxel grid for base volume (0-based indexing)
    [Yi, Xi, Zi] = ndgrid(0:size(baseVol,1)-1, ...
                          0:size(baseVol,2)-1, ...
                          0:size(baseVol,3)-1);

    baseVox = [Yi(:), Xi(:), Zi(:), ones(numel(Yi), 1)];

    % Convert base voxel indices to world coordinates
    worldPts = baseVox * A_base';

    % Convert world coordinates to moving image voxel indices
    movVox = worldPts * inv(A_mov)';

    ym = movVox(:,1) + 1;
    xm = movVox(:,2) + 1;
    zm = movVox(:,3) + 1;

    
    % Interpolate moving volume at computed voxel positions
    volOnBase = interpn(double(movVol), ym, xm, zm, method, NaN);
    % spline kernels smear a few epsilon-level numbers into the NaN area
    volOnBase( abs(volOnBase) < 1e-5 ) = NaN;
    volOnBase = reshape(volOnBase, size(baseVol));
end

function A = getAffineMatrix(info)
% Constructs a 4x4 voxel-to-world affine matrix from NIfTI header information

    % Check for sform matrix
    if isfield(info, 'raw') && isfield(info.raw, 'sform_code') && info.raw.sform_code > 0
        % Construct affine matrix from srow_x, srow_y, srow_z
        A = [info.raw.srow_x; info.raw.srow_y; info.raw.srow_z; 0 0 0 1];
    elseif isfield(info, 'raw') && isfield(info.raw, 'qform_code') && info.raw.qform_code > 0
        % Construct affine matrix from qform parameters
        b = info.raw.quatern_b;
        c = info.raw.quatern_c;
        d = info.raw.quatern_d;
        a = sqrt(1.0 - (b^2 + c^2 + d^2));

        qfac = info.raw.pixdim(1);
        if qfac == 0
            qfac = 1;
        end

        R = [a^2 + b^2 - c^2 - d^2,     2*b*c - 2*a*d,         2*b*d + 2*a*c;
             2*b*c + 2*a*d,             a^2 + c^2 - b^2 - d^2, 2*c*d - 2*a*b;
             2*b*d - 2*a*c,             2*c*d + 2*a*b,         a^2 + d^2 - b^2 - c^2];

        % Apply voxel dimensions
        pixdim = info.raw.pixdim(2:4);
        R = R * diag(pixdim);

        % Apply qfac to the third column
        R(:,3) = R(:,3) * qfac;

        % Construct affine matrix
        A = eye(4);
        A(1:3,1:3) = R;
        A(1:3,4) = [info.raw.qoffset_x; info.raw.qoffset_y; info.raw.qoffset_z];
    else
        % Fallback: use Transform.T if available
        if isfield(info, 'Transform') && isfield(info.Transform, 'T')
            A = info.Transform.T;
        else
            error('No valid affine transformation found in NIfTI header.');
        end
    end
end
