function spmup_plotmotion(file)

% simple routine to plot motion parameters in mm and degrees
% also added is the total displacement
%
% FORMAT: spmup_plotmotion(file)
%
% INPUT: the file in which the motion paramters are
%
% OUTPUT: a figure 
%
% ------------------------------
% Copyright (C) spmup team 2019

%% get the data
current = pwd;
if nargin == 0
    motion_file = uigetfile(pwd,'select realignemt parameter file');
    if filepath == 0
        return
    else
        cd(filepath);
    end

end

%% compute
motion_param       = load(motion_file.name);
[n,~]              = size(motion_param);
temp_motion        = motion_param; 
temp_motion(:,4:6) = temp_motion(:,4:6)*180/pi; % ~57mm head size
derivatives        = [0 0 0 0 0 0 ; diff(temp_motion)];
RMS                = sqrt(sum(derivatives.^2,2);  % Mean square displacement 

%% plot
figure('Name','Motion','Visible','On');
set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized', 'outerposition',[0 0 1 1]);
subplot(2,2,1); plot(temp_motion(:,1:3),'LineWidth',3);
axis tight; grid on; title('Translation','Fontsize',14)
xlabel('scans'); ylabel('mm');
subplot(2,2,2); plot(temp_motion(:,4:6),'LineWidth',3);
axis tight; grid on; title('Rotation','Fontsize',14)
xlabel('scans'); ylabel('degrees');
subplot(2,2,3:4); plot(delta,'LineWidth',3); hold on;
axis tight; grid on; title(['Mean square displacement ' num2str(mean(RMS)) ' std ' num2str(std(RMS))],'FontSize',14);
xlabel('scans'); ylabel('displacement');
saveas(gcf, 'motion plots.fig','fig');
print (gcf,'-dpsc2', '-bestfit', [pwd filesep 'motion plots.ps']);
close(gcf)
cd(current)
