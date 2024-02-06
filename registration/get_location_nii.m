function location_info = get_location_nii(niiHdr)
dims = niiHdr.dime;
hst = niiHdr.hist;
R = get_rotation_nifti(niiHdr);
R(1,1)=-R(1,1);
R(2,2)=-R(2,2);
R(3,3)=dims.pixdim(1)*R(3,3);

location_info = cell(5,2);
location_info(:,1)={'Size';'Spacing';'Origin';'Index';'Direction'};

location_info{1,2}= num2str(dims.dim(2:4));
location_info{2,2}= num2str(dims.pixdim(2:4));
location_info{3,2}= num2str([-hst.qoffset_x -hst.qoffset_y hst.qoffset_z]);

location_info{4,2}= num2str([0 0 0]);
location_info{5,2}= num2str(R(:)');
end