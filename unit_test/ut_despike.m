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
            FilteredM(index,:) = medfilt1(squeeze(Data(x,y,z,:)),[],'truncate',3); 
        end
    end
end

flags = struct('auto_mask','off','method','median','window',3);
[despiked,filtered]  = spmup_despike(Data,[],flags);
index = 1;
for x=1:5
    for y=1:5
        for z=1:5
            despikedM(index,:) = squeeze(despiked(x,y,y,:)); 
            filteredM(index,:) = squeeze(filtered(x,y,y,:)); 
            index = index+1;
        end
    end
end

%% do some plots
figure
subplot(2,2,1);
plot(squeeze(DataM(1,:)),'*b','LineWidth',2);
hold on; grid on; axis tight
plot(squeeze(filteredM(1,:)),'r','LineWidth',2);
title('Exemple of filter with spmup despike vs data');

[~,position] = max(sqrt(mean([DataM-filteredM].^2,2)));
subplot(2,2,2);
plot(DataM(position,:),'*','LineWidth',2);
hold on; grid on; axis tight
plot(FilteredM(position,:),'k','LineWidth',2);
plot(filteredM(position,:),'--r','LineWidth',2);
title('Most diff data to filtered version (medfilt1 / spmup)');

subplot(2,2,3);
plot(median(DataM,1),'*','LineWidth',2);
hold on; grid on; axis tight
plot(median(despikedM,1),'r','LineWidth',2);
title('mean time course of desipked data');
subplot(2,2,4);
plot(median(DataM,1),'*','LineWidth',2);
hold on; grid on; axis tight
plot(nanmedian(FilteredM,1),'k','LineWidth',3);
plot(nanmedian(filteredM,1),'--r','LineWidth',2);
title('Data with filtered versions (medfilt1 / spmup)');

%% quantify
D = (FilteredM == filteredM);
fprintf('%g%% time courses identical to median filter (except end)\n',mean(mean(D(:,1:100)<eps))*100);
M(:,1) = sqrt(mean([DataM(indices,:) - FilteredM(indices,:)].^2,2));
M(:,2) = sqrt(mean([DataM(indices,:) - filteredM(indices,:)].^2,2));

overfit = sum(M(:,1) == 0) / size(M,1);
fprintf('matlab median filter matches exactly the data in %g%% of cases \n',overfit*100)
goodfit = sum(M(:,2) < 1) / size(M,1);
fprintf('spmup filter matches the data correctly in %g%% of cases (rms<1)\n',goodfit*100)
bad_fit = find(M(:,2)>=1);
fprintf('spmup filter gives a bad fit in %g%% of all data \n',(numel(Data)-length(bad_fit))/numel(Data));

figure; 
for i=1:length(bad_fit)
    subplot(1,2,1);
    plot([1:101],DataM(bad_fit(i),:),'*b','LineWidth',2);
    hold on; grid on; axis tight
    plot([1:101],DataM(bad_fit(i),:),'b','LineWidth',2);
    plot(despikedM(bad_fit(i),:),'r','LineWidth',2);
    title(['filtered data ' num2str(indices(i)) ' RMS ' num2str(M(bad_fit(i),2))]);
    hold off

    subplot(1,2,2);
    plot([1:101],DataM(bad_fit(i),:),'*b','LineWidth',2);
    hold on; grid on; axis tight
    plot([1:101],DataM(bad_fit(i),:),'b','LineWidth',2);
    plot(despikedM(bad_fit(i),:),'r','LineWidth',2);
    title(['despiked data ' num2str(indices(i)) ' RMS ' num2str(sqrt(mean([DataM(bad_fit(i),:) - despikedM(bad_fit(i),:)].^2,2)))]);
    pause; hold off
end

%% do additional testing for other smoother - just to check the function works
 flags.method = 'moving';
[despiked2,filtered2]  = spmup_despike(Data,[],flags);

