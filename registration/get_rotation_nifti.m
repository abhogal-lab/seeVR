function R = get_rotation_nifti(niiHdr)

b = niiHdr.hist.quatern_b;
c = niiHdr.hist.quatern_c;
d = niiHdr.hist.quatern_d;

% Sometimes rounding results in negative values for 1-b^2-c^2-d^2,
% therefore the max is in there.
a = sqrt(max(1-b^2-c^2-d^2,0));

R =    [a*a+b*b-c*c-d*d   2*b*c-2*a*d       2*b*d+2*a*c     ;...
        2*b*c+2*a*d       a*a+c*c-b*b-d*d   2*c*d-2*a*b     ;...
        2*b*d-2*a*c       2*c*d+2*a*b       a*a+d*d-c*c-b*b ];