%% unit testing for spmup_despike
% uses medfilt1(x,3) from matlab on a small synthetic 4D set 
% the comparison focuses on the median filter itself from which one derives
% the despiked data - but if the filter is bad so will be the despiked data
% outputs the percentage difference, look at what is different and by how
% much -- on average ~80% are identifcal and differing are <1% of the
% initial data range ; I think that good enough since the differences are
% also symetric - once averaged in time it's almost null
% ------------------------------
% Cyril Pernet 17 Decembre 2016

%% make data up

DataM = NaN(125,101);
Data = NaN(5,5,5,101);
FilteredM = NaN(125,101);
FilteredData = NaN(5,5,5,101);

fs = 100;     % Sampling rate
t = 0:1/fs:1; % Time vector
% Noise Signal - Input
index = 1;
for x=1:5
    for y=1:5
        for z=1:5
            DataM(index,:) = sin(2*pi*t*randi(10,1,1))+.25*sin(2*pi*t*randi(50,1,1)); 
            Data(x,y,z,:) =  DataM(index,:); 
            FilteredM(index,:) = medfilt1(squeeze(Data(x,y,z,:)),3); % Median filtering - Output
            FilteredData(x,y,z,:) = FilteredM(index,:) ; index = index+1;
        end
    end
end

flags = struct('auto_mask','off','method','median','window',3);
[despiked,filtered]  = spmup_despike(Data,[],flags);

%% do some plots
figure
subplot(2,2,1);
plot(squeeze(Data(1,1,1,:)),'*','LineWidth',2);
hold on; grid on; axis tight
plot(squeeze(despiked(1,1,1,:)),'r','LineWidth',2);
title('1 time course with spmup desipked version');
subplot(2,2,2);
plot(squeeze(Data(1,1,1,:)),'*','LineWidth',2);
hold on; grid on; axis tight
plot(squeeze(FilteredData(1,1,1,:)),'k','LineWidth',3);
plot(squeeze(filtered(1,1,1,:)),'--r','LineWidth',2);
title('Data with filtered versions (medfilt1 / spmup)');

index = 1;
for x=1:5
    for y=1:5
        for z=1:5
            despikedM(index,:) = squeeze(despiked(x,y,y,:)); 
            filteredM(index,:) =  squeeze(filtered(x,y,y,:)); 
            index = index+1;
        end
    end
end

subplot(2,2,3);
plot(mean(DataM,1),'*','LineWidth',2);
hold on; grid on; axis tight
plot(mean(despikedM,1),'r','LineWidth',2);
title('mean time course with spmup desipked version');
subplot(2,2,4);
plot(mean(DataM,1),'*','LineWidth',2);
hold on; grid on; axis tight
plot(nanmean(FilteredM,1),'k','LineWidth',3);
plot(nanmean(filteredM,1),'--r','LineWidth',2);
title('Data with filtered versions (medfilt1 / spmup)');

%% quantify
D = (filteredM == FilteredM);
fprintf('%g time courses with the same time course (except start)\n',sum(sum(D,2)<=2)/125*100);
figure; indices = find(sum(D,2)>2);
for i=1:length(indices)
    plot([1:101],DataM(indices(i),:),'*','LineWidth',2);
    hold on; grid on; axis tight
    plot(FilteredM(indices(i),:),'k','LineWidth',3);
    plot(filteredM(indices(i),:),'--r','LineWidth',2);
    title(['filtered data ' num2str(indices(i))]);
    pause; hold off
end

figure; 
D = (filteredM - FilteredM);
D = mean(D(indices,:),1) / mean(range(DataM(indices,:),2)) .* 100;
plot(D,'LineWidth',3);
title('average difference spmup - medfilt1 among different time courses')
hold on; grid on; axis tight; ylabel('difference in filter (unit % data range)')
