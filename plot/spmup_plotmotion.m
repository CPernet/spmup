function spmup_plotmotion(folder)

% simple routine to plot motion parameters in mm and degrees
% also added is the total displacement
%
% FORMAT: spmup_plotmotion(folder)
%
% INPUT: the folder in which the motion paramters are
%        if not providing, search the current folder
%
% OUTPUT: a figure (both .fig and .jpg)
%
% Cyril Pernet - 02 Dec 2016
% ------------------------------
% Copyright (C) spmup team 2016

%% get the data
current = pwd;
if nargin == 0
    filepath = uigetdir(pwd,'select fmri run directory');
    if filepath == 0
        return
    else
        cd(filepath);
    end
elseif nargin == 1
    cd(varargin{1});
end

motion_file = dir('rp*.txt'); % SPM
if isempty(motion_file)
    motion_file = dir('*_mcf.par'); % FSL
end

if isempty(motion_file)
    cd(current); error('no motion paramter files found')
end

%% compute
motion_param       = load(motion_file.name);
[n,~]              = size(motion_param);
temp_motion        = motion_param; 
temp_motion(:,4:6) = temp_motion(:,4:6)*180/pi; % use deg rather than radian
derivatives        = diff(temp_motion);

delta = zeros(n,1);  % Mean square displacement in two scans
for i = 2:n
    delta(i) = sqrt((temp_motion(i-1,1) - temp_motion(i,1))^2 + ...
        (temp_motion(i-1,2) - temp_motion(i,2))^2 +...
        (temp_motion(i-1,3) - temp_motion(i,3))^2 +...
        (temp_motion(i-1,4) - temp_motion(i,4))^2 +...
        (temp_motion(i-1,5) - temp_motion(i,5))^2 +...
        (temp_motion(i-1,6) - temp_motion(i,6))^2);
end

%% plot
figure('Name','Motion outlier detection','Visible','On');
set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized', 'outerposition',[0 0 1 1]);
subplot(2,2,1); plot(temp_motion(:,1:3),'LineWidth',3);
axis tight; grid on; title('Translation','Fontsize',14)
xlabel('scans'); ylabel('mm');
subplot(2,2,2); plot(temp_motion(:,4:6),'LineWidth',3);
axis tight; grid on; title('Rotation','Fontsize',14)
xlabel('scans'); ylabel('degrees');
subplot(2,2,3:4); plot(delta,'LineWidth',3); hold on;
axis tight; grid on; title(['Mean square displacement ' num2str(mean(delta)) ' std ' num2str(std(delta))],'FontSize',14);
xlabel('scans'); ylabel('displacement');
saveas(gcf, 'motion plots.fig','fig');
print (gcf,'-dpsc2', '-bestfit', [pwd filesep 'motion plots.ps']);
close(gcf)

cd(current)
