function spmup_hrf_boost(varargin)

% this function takes a SPM.mat and related files
% to create boosted beta parameters and con images
% if the hrf doesn't fit properly the data, some
% information can be recovered using the 1st and 2nd 
% derivatives - 
%
% FORMAT spm_hrf_boost
%        spm_hrf_boost(SPM_location,shift)
%
% INPUT SPM_location is the full name of the SPM.mat
%       shift is time shift allowed around the hrf peak - default is 2.5
%
% OUTPUT created in /hrf_boost
%        beta_XXXX corresponding to the boosted hrf
%        con_XXXX combined boosted hrf regressors
%
% Cyril Pernet April 2014
% Updated August 2015 to handle .nii + batch input
% -----------------------------------------------
% Copyright (C) spmup team 2015

defaults = spm_get_defaults;
img_ext = defaults.images.format;

current = pwd;
shift = 2;
if nargin == 0
    [t,sts] = spm_select(1,'mat','Select 1st level SPM file');
    if sts == 0
        return
    end
elseif nargin == 1
    t = varargin{1};
elseif nargin == 2
    t = varargin{1};
    shift = varargin{2};
end

if iscell(t); t = cell2mat(t); end
if iscell(shift); shift = cell2mat(t); end
if isstruct(shift); shift = getfield(shift,cell2mat(fieldnames(shift))); end

% check the design 
% ------------------
load(t);

if strcmp(SPM.xBF.name,'hrf (with time and dispersion derivatives)')
    type = 3;
elseif strcmp(SPM.xBF.name,'hrf (with time derivative)')  
    type = 2;
else
    error('spm_hrf_boost only works for designs with 1st and/or 2nd derivatives')
end

[path,file,ext] = fileparts(t);
cd(path);
mkdir('hrf_boost')

% using a reference trial, estimate the range of possible delays
% use this to constrains the hrf boost, since we don't want to boost
% the hrf with a time derivative that fits some artefacts - 
% take as possible range +/- the shift parameter
% --------------------------------------------------------------
disp('--------------------------------')
disp('  Estimating acceptable delays  ')
disp('--------------------------------')

% parameters of the design
xBF.dt = SPM.xBF.dt;
xBF.name = SPM.xBF.name;
xBF.length = SPM.xBF.length;
xBF.order = SPM.xBF.order;
xBF = spm_get_bf(xBF);

% orthogonalise derivative(s) and normalize
xBF.bf  =  spm_orth(xBF.bf);
SS = xBF.bf'*xBF.bf;
for f=1:size(xBF.bf,2)
    xBF.bf(:,f) = xBF.bf(:,f)./sqrt(SS(f,f));
end

% the range of time covered by this model is
BFtime = [0:xBF.dt:(length(xBF.bf)-1)*xBF.dt];
PeakTime = BFtime(find(xBF.bf(:,1) == max(xBF.bf(:,1))));
PeakDelays = [PeakTime-shift PeakTime+shift];
fprintf('parameters will be boosted for estimated responses peaking between %g and %g sec \n',PeakDelays(1),PeakDelays(2));

% boost parameters 
% ------------------

disp('--------------------------------')
disp('  Re-estimating hrf magnitude   ')
disp('--------------------------------')

index = 1; i = 1;
for s=1:size(SPM.Sess,2) % for each session
    for n = 1:length(SPM.Sess(s).U) % for each condition
        for c = 1:length(SPM.Sess(s).U(n).name) % for each regressor in this condition
            
            % name
            if index < 10
                name = sprintf('boost_beta_000%g',index);
                im(1,:) = [pwd filesep sprintf('beta_000%g.%s',index,img_ext)];
                im(2,:) = [pwd filesep sprintf('beta_000%g.%s',index+1,img_ext)];
                if type == 3
                    im(3,:) = [pwd filesep sprintf('beta_000%g.%s',index+2,img_ext)];
                end
            elseif (index>=10) && (index<100)
                name = sprintf('boost_beta_00%g',index,img_ext);
                im(1,:) = [pwd filesep sprintf('beta_00%g.%s',index,img_ext)];
                im(2,:) = [pwd filesep sprintf('beta_00%g.%s',index+1,img_ext)];
                if type == 3
                    im(3,:) = [pwd filesep sprintf('beta_00%g.%s',index+2,img_ext)];
                end
            elseif (index>=100) && (index<1000)
                name = sprintf('boost_beta_0%g',index,img_ext);
                im(1,:) = [pwd filesep sprintf('beta_0%g.%s',index,img_ext)];
                im(2,:) = [pwd filesep sprintf('beta_0%g.%s',index+1,img_ext)];
                if type == 3
                    im(3,:) = [pwd filesep sprintf('beta_0%g.%s',index+2,img_ext)];
                end
            elseif (index>=1000) && (index<10000)
                name = sprintf('boost_beta_%g',index,img_ext);
                im(1,:) = [pwd filesep sprintf('beta_%g.%s',index,img_ext)];
                im(2,:) = [pwd filesep sprintf('beta_%g.%s',index+1,img_ext)];
                if type == 3
                    im(3,:) = [pwd filesep sprintf('beta_%g.%s',index+2,img_ext)];
                end
            end
            
            % dirty update for these special cases
            if index == 9 
                im(1,:) = [pwd filesep sprintf('beta_000%g.%s',index,img_ext)];
                im(2,:) = [pwd filesep sprintf('beta_00%g.%s',index+1,img_ext)];
                if type == 3
                    im(3,:) = [pwd filesep sprintf('beta_00%g.%s',index+2,img_ext)];
                end
            end
 
            if index == 99 
                im(1,:) = [pwd filesep sprintf('beta_00%g.%s',index,img_ext)];
                im(2,:) = [pwd filesep sprintf('beta_0%g.%s',index+1,img_ext)];
                if type == 3
                    im(3,:) = [pwd filesep sprintf('beta_0%g.%s',index+2,img_ext)];
                end
            end
            
            if index == 999
                im(1,:) = [pwd filesep sprintf('beta_0%g.%s',index,img_ext)];
                im(2,:) = [pwd filesep sprintf('beta_%g.%s',index+1,img_ext)];
                if type == 3
                    im(3,:) = [pwd filesep sprintf('beta_%g.%s',index+2,img_ext)];
                end
            end
            
            % inform user
            if type == 2
                fprintf('Combining \n %s \n %s \n',im(1,:),im(2,:));
            elseif type == 3
                fprintf('Combining \n %s \n %s \n %s \n',im(1,:),im(2,:),im(3,:));
            end
            
            % Create a mask image for voxels outside the PeakDelays
            V = spm_vol(im);
            sts = spm_check_orientations(V);
            images = spm_read_vols(V);
            T2P = NaN(V(1).dim); % time to peak image
            voxels = find(~isnan(squeeze(images(:,:,:,1))));
            f = waitbar(0,'Estimating time to peak voxel-wise','name',['HRF BOOST ' num2str(index)]);
            for v=1:length(voxels)
                waitbar(v/length(voxels));
                [x,y,z] = ind2sub(size(squeeze(images(:,:,:,1))),voxels(v));
                
                if type == 2
                    Model = SPM.xBF.bf * [images(x,y,z,1) images(x,y,z,2)]';
                elseif type == 3
                    Model = SPM.xBF.bf * [images(x,y,z,1) images(x,y,z,2) images(x,y,z,3)]';
                end
                
                if images(x,y,z,1) > 0
                    T2P(voxels(v)) = BFtime(find(Model==max(Model)));
                else
                    T2P(voxels(v)) = BFtime(find(Model==min(Model)));
                end
            end
            close(f)
            newV = spm_create_vol(V(1));
            [path,file,ext] = fileparts(V(1).fname);
            newV.fname = [path filesep 'hrf_boost' filesep 'T2P' num2str(index) ext];
            newV.descrip = sprintf('time to peak image of %s',V(1).descrip);
            out = spm_write_vol(newV,T2P);
            Mask = logical((T2P>=PeakDelays(1)).*(T2P<=PeakDelays(2)));
%            Mask = ~isnan(squeeze(images(:,:,:,1)));
%             figure; for z=1:V(1).dim(3)
%                 imagesc(squeeze(Mask(:,:,z))); pause(0.2);
%             end
            
            % boost the hrf values
            X = SPM.xX.X(:,[index index+1 index+2]);
            SumOfSquares = diag(X'*X); % these are the weights
            if sum(SumOfSquares - 1) > 0.0001 % also check X is normalized
                X(:,1) = X(:,1)./norm(X(:,1));
                X(:,2) = X(:,2)./norm(X(:,2));
                X(:,3) = X(:,3)./norm(X(:,3));
                SumOfSquares = diag(X'*X);
            end
            
            % Create and execute the spm image calculation.
            if type == 2
                magnitude = sqrt(((squeeze(images(:,:,:,1)).^2).*SumOfSquares(1)) + ...
                ((squeeze(images(:,:,:,2)).^2).*SumOfSquares(2))); % like Steffener
            elseif type == 3
                magnitude = sqrt(((squeeze(images(:,:,:,1)).^2).*SumOfSquares(1)) + ...
                ((squeeze(images(:,:,:,2)).^2).*SumOfSquares(2)) + ...
                ((squeeze(images(:,:,:,3)).^2).*SumOfSquares(3))); 
            end
            sign = (squeeze(images(:,:,:,1))) ./ abs((squeeze(images(:,:,:,1)))); % get sign like Calhoun
            boosted_values = magnitude .* sign;
            
            % Only change values within Mask
            boosted_hrf = squeeze(images(:,:,:,1)); % beta hrf
            boosted_hrf(Mask) = boosted_values(Mask);
%             figure; for z=1:V(1).dim(3)
%                 imagesc(boosted_hrf(:,:,z)); pause(0.2);
%             end
            
            % save file
            newV = spm_create_vol(V(1));
            [path,file,ext] = fileparts(V(1).fname);
            newV.fname = [path filesep 'hrf_boost' filesep name ext];
            newV.descrip = sprintf('Boosted version of %s',V(1).descrip);
            out = spm_write_vol(newV,boosted_hrf);
            fprintf('boost_%s done \n',file);
            
            % update index for the next hrf
            hrf_indices(i) = index; i=i+1;
            if type == 2
                index = index+2;
            elseif type == 3
                index = index+3;
            end
            clear name im V images X magnitude sign boosted_hrf newV
        end
    end
    % after a given session, we might have to move further if motion etc
    % are added to the design
    index = index + size(SPM.Sess(s).C.C,2);
end

% also update constrasts
% ----------------------

if isfield(SPM,'xCon')
    
    disp(' ')
    disp('--------------------------------')
    disp('    Re-estimating contrasts     ')
    disp('--------------------------------')
    cd hrf_boost/
    
    for c = 1:length(SPM.xCon) % for each contrast
        columns = find(SPM.xCon(c).c); % check columns
        test = intersect(columns,hrf_indices);
        if length(test) == length(columns) % constrast involving combination of hrf only
            % load boosted parameters
            for i=1:length(columns)
                if columns(i) < 10
                    im(i,:) = [pwd filesep sprintf('boost_beta_000%g.%s',columns(i),img_ext)];
                elseif columns(i) < 10
                    im(i,:) = [pwd filesep sprintf('boost_beta_00%g.%s',columns(i),img_ext)];
                elseif columns(i) < 100
                    im(i,:) = [pwd filesep sprintf('boost_beta_0%g.%s',columns(i),img_ext)];
                end
            end
            V = spm_vol(im);
            sts = spm_check_orientations(V);
            images = spm_read_vols(V);
            
            % compute the constrast
            fprintf('Contrasting \n')
            boosted_con = zeros(size(images,1),size(images,2),size(images,3));
            for i=1:length(columns)
                fprintf('%s \n',im(i,:));
                boosted_con = boosted_con + (SPM.xCon(c).c(columns(i)).*squeeze(images(:,:,:,i)));
            end
            
            % save file
            newV = spm_create_vol(V(1));
            [path,file,ext] = fileparts(SPM.xCon(c).Vcon.fname);
            newV.fname = [pwd filesep 'boost_' file ext];
            newV.descrip = sprintf('Boosted version of %s',SPM.xCon(c).Vcon.descrip);
            out = spm_write_vol(newV,boosted_con);
            fprintf('boost_%s done \n',SPM.xCon(c).Vcon.fname);
        end
    end
end

cd(current)
disp('------------------')
disp('spm_hrf_boost done')

        
        
        
        