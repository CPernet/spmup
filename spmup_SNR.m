function spmup_SNR(image)

% simple routine to compute the SNR from a single 3D image
% inspired largley by nii_mean_stdev Chris Rorden, 2014
% http://opensource.org/licenses/BSD-2-Clause

[hdr, meanImg, sdImg, n] = statSub3D(normBrightness, fnms);
[hdr, meanImg, sdImg, n] = sumSub(hdr, meanImg, sdImg, n, normBrightness, fnms);
sdImg = sqrt(sdImg/(n-1)); %convert to standard deviation
for i=1:size(fnms,1)
	fnm = deblank(fnms(i,:));
    [hdr, img] = read_volsSub (fnm);
    if normBrightness, img = normBrightnessSub(img); end;
    if isempty(meanImg), meanImg = zeros(size(img)); end;
    if isempty(m2Img), m2Img = zeros(size(img));  end;
    n = n + 1;
    delta = img - meanImg;
    meanImg = meanImg + (delta / n);
    m2Img = m2Img + delta.*(img-meanImg);
    %fprintf('%f\t%f\t%f\n',img(1),meanImg(1),m2Img(1));
end
function img = normBrightnessSub(img)
img = img - min(img(:));
mdn = median(img(img > 0));
if mdn == 0, return; end;
img = img / mdn;

