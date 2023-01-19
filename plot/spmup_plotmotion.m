function spmup_plotmotion(motion_file,radius)

% simple routine to plot motion parameters in mm and degrees
% also added is the root mean square displacement
%
% FORMAT: spmup_plotmotion(motion_param,radius)
%
% INPUT: motion_param is the file in which the motion paramters are
%        radius the head size (useful to compute displacement) if empty use 57 mm
%
% OUTPUT: a figure 
%
% Cyril Pernet
% --------------------------
%  Copyright (C) SPMUP Team 

%% get the data
if nargin == 0
    motion_file = uigetfile(pwd,'select realignemt parameter file');
    if filepath == 0
        return
    else
        cd(filepath);
    end
end

if ~exist('radius','var')
    radius = 180/pi; % ~57mm head size
end

%% compute
motion_param       = load(motion_file.name);
temp_motion        = motion_param; 
temp_motion(:,4:6) = temp_motion(:,4:6)*radius; 
derivatives        = [0 0 0 0 0 0 ; diff(temp_motion)];
RMS                = sqrt(sum(derivatives.^2,2));  % Mean square displacement 

%% plot
figure('Name','Motion','Visible','On');
set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized', 'outerposition',[0 0 1 1]);

subplot(2,2,1); plot(temp_motion(:,1:3),'LineWidth',3);
axis tight; grid on; title('Translation','Fontsize',14)
xlabel('scans'); ylabel('mm');

subplot(2,2,2); plot(temp_motion(:,4:6),'LineWidth',3);
axis tight; grid on; title('Rotation','Fontsize',14)
xlabel('scans'); ylabel('degrees');

subplot(2,2,3:4); plot(RMS,'LineWidth',3); hold on;
axis tight; grid on; title(['Mean square displacement ' num2str(mean(RMS)) ' std ' num2str(std(RMS))],'FontSize',14);
xlabel('scans'); ylabel('displacement');

