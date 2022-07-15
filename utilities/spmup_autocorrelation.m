function R = spmup_autocorrelation(varargin)

% this function computes efficiently the autocorrelation of all in mask 
% voxels of time series (fMRI) data.
%
% FORMAT R = spmup(timeseries,mask) 
%
% INPUT timeseries can be:
%       - a structure array containing image volume information (see spm_vol)
%       - a 2D matrix with observation in rows and time in columns
%       mask (optional)
%       - a struture for an image maching timeseries, allows restricting
%       analyses to some voxels (by default all non 0 time courses are analzed)
%
% OUTPUT R is a 2D or 3D matrix of autocorrelation coefficients
%
% The Wiener-Khinchin theorem relates the autocorrelation function to the 
% power spectral density via the Fourier transform. It states that the 
% autocorrelation function of a wide-sense-stationary random process 
% (the 1st moment does not vary in time) has a spectral decomposition given 
% by the power spectrum of that process. Since FFT can be applied to a
% matrix, this means we can get all autocorrelation functions estimates 
% in one pass - and then derive the coefficents, i.e. the distance in the 
% spectrum between the peak and it's neighbourgh (harmonic?) which is 
% equivalent to the lag of the autocorrelation.
%
% Cyril Pernet v1: 13 Decembre 2016
% --------------------------------------------------------------------------
% Copyright (c) SPM Utility Plus toolbox

%% deal wit the data in

if nargin == 2
    try 
        if ischar(varargin{2})
            index = find(spm_read_vols(varargin{2})); % in mask voxel
        else
            index = find(varargin{2});
        end
        [X, Y, Z] = ind2sub(varargin{1}(1).dim,index); % coordinates
        Data      = spm_get_data(varargin{1},[X Y Z]'); % 2D matrix
        Format    = 'image';
    catch
        error('could not read the data')
    end
else
    Data   = varargin{1};
    Format = 'matrix';
    clear varargin
end

% extra check
if length(size(Data)) ~= 2
    error('the data used don''t seem to make a 2D array');
else
    [p,n] = size(Data);
    if isvector(Data) && n>1
        Data  = Data';
        [p,n] = size(Data);
    end
end


%% compute the autocorrelation using fft
disp('estimating autocorrelation for each voxel')

nfft = 2^nextpow2(2*p-1);
Data = spm_detrend(Data,1);
R    = ifft(fft(Data,nfft).*conj(fft(Data,nfft)));
R    = R(1:p,:)'; % [R(end-p+2:end,:) ; R(1:p,:)]'; for a full window

%% now the slow part
autocorrelation_window = NaN(n,1);

for v=1:n
    if sum(R(v,:)) == 0
        autocorrelation_window(v) = 0;
    else
        [~,locs] = findpeaks(R(v,:)); % get local maxima
        if ~isempty(locs)
            [~,maxloc]                = max(R(v,:));   % find spectral peak
            tmp                       = locs-maxloc;
            tmp(tmp <= 0)             = NaN;           % avoid being on itself
            [~,closestpeak]           = min(tmp);      % find next peak
            closestpeak               = locs(closestpeak);
            [~,minloc]                = min(R(v,maxloc:closestpeak));
            autocorrelation_window(v) = maxloc+minloc; % return the distance
        else
            autocorrelation_window(v) = 0;
        end
    end
end

%% reformat as input
clear R
if strcmp(Format,'image')
    R        = NaN(varargin{1}(1).dim);
    R(index) = autocorrelation_window;
else
    R        = autocorrelation_window;
end

