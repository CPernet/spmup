function MakeContours(xSPM)

% simple routine to save the contours of activated clusters as a binary
% image from the current plotted data (i.e. using xSPM)
%
% the routine uses a simple cube as structuring element for the dilation
% so the output is roughly of the specified contour size
%
% Cyril Pernet 26 Dec 2013

if nargin == 0
    try
        xSPM = evalin('base', 'xSPM');
    catch NARGIN_ERROR
        error('Please run SPM results query first')
    end
end

[contour_size] = spm_input('contours size (in voxel)',1,'i',' ',1);



%-Get filename
%--------------------------------------------------------------------------
F   = xSPM.title;
F   = [F '_contour.img'];
spm('Pointer','Watch')

%-Set up header information
%--------------------------------------------------------------------------
Vo  = struct(...
        'fname',    F,...
        'dim',      xSPM.DIM',...
        'dt',       [spm_type('float32') spm_platform('bigend')],...
        'mat',      xSPM.M,...
        'descrip',  'binary contour map');
Vo.dt(1) = spm_type('uint8');

%-Reconstruct (filtered) image from XYZ & Z pointlist
%--------------------------------------------------------------------------
Y      = nan(xSPM.DIM(1:3)');
OFF    = xSPM.XYZ(1,:) + xSPM.DIM(1)*(xSPM.XYZ(2,:)-1 + xSPM.DIM(2)*(xSPM.XYZ(3,:)-1));
Y(OFF) = xSPM.Z > 0;
Y(find(isnan(Y(:)))) = 0; % Y is a binary of the clusters
% figure; for z=1:xSPM.DIM(3); imagesc(Y(:,:,z)); pause(0.2); end

% structuring element of dilation
square = ones(2+contour_size,2+contour_size);
se = strel('arbitrary',square,square.*(2+contour_size));

% apply
Y2 = imdilate(Y,se);
% figure; for z=1:xSPM.DIM(3); imagesc(Y2(:,:,z)); pause(0.2); end

% substract to get contour
Y3 = Y2-Y; Y3 = Y3-(min(Y3(:)));
for z=1:xSPM.DIM(3)
    index(z) = range(range(Y3(:,:,z)));
end

figure; for z=find(index); imagesc(Y3(:,:,z)); pause(0.2); end
close


%-Write the reconstructed volume
%--------------------------------------------------------------------------
Vo = spm_write_vol(Vo,Y3);
spm('alert"',{'Written:',['    ',spm_select('CPath',F)]},mfilename,1);

