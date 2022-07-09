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

clear variables
Data     = NaN(5,5,5,101);
Filtered = Data;
fs       = 100;     % Sampling rate
t        = 0:1/fs:1; % Time vector

% Noise Signal - Input
for x=1:5
    for y=1:5
        for z=1:5
            Data(x,y,z,:)     = sin(2*pi*t*randi(10,1,1))+.25*sin(2*pi*t*randi(50,1,1));
            Filtered(x,y,z,:) = medfilt1(squeeze(Data(x,y,z,:)),[],'truncate',3);
        end
    end
end

flags = struct('auto_mask','off','method','median','window',3,'savelog','off');
[despiked,filtered]  = spmup_despike(Data,[],flags);

index = 1;
for x=1:5
    for y=1:5
        for z=1:5
            filteredM(index,:) = squeeze(filtered(x,y,z,:));
            FilteredM(index,:) = squeeze(Filtered(x,y,z,:));
            DM(index,:) = [Data(x,y,z,:)-Filtered(x,y,z,:)].^2; % difference to matlab
            DS(index,:) = [Data(x,y,z,:)-filtered(x,y,z,:)].^2; % difference to SPMUP
            D(index,:)  = squeeze(Data(x,y,z,:));
            S(index,:)  = squeeze(despiked(x,y,z,:));
            index = index+1;
        end
    end
end

%% do some plots
figure
[~,position] = min(sqrt(mean(DS,2)));
subplot(2,2,1);
plot(squeeze(D(position,:)),'*b','LineWidth',2);
hold on; grid on; axis tight
plot(squeeze(filteredM(position,:)),'-r','LineWidth',2);
title('Best fit spmup despike vs data');

[~,position] = max(sqrt(mean(DS,2)));
subplot(2,2,2);
plot(D(position,:),'*','LineWidth',2);
hold on; grid on; axis tight
plot(FilteredM(position,:),'k','LineWidth',2);
plot(filteredM(position,:),'--r','LineWidth',2);
title('Most diff data to filtered version (medfilt1 / spmup)');

subplot(2,2,3);
plot(mean(D,1),'*','LineWidth',2);
hold on; grid on; axis tight
plot(mean(S,1),'r','LineWidth',2);
title('mean time course of desipked data');

subplot(2,2,4);
plot(mean(D,1),'*','LineWidth',2);
hold on; grid on; axis tight
plot(nanmean(FilteredM,1),'k','LineWidth',3);
plot(nanmean(filteredM,1),'--r','LineWidth',2);
title('Data with filtered versions (medfilt1 / spmup)');

%% quantify
DF = (FilteredM(:,1:100) == filteredM(:,1:100));
fprintf('%g%% time courses identical to median filter \n',mean(mean(DF,2)<eps)*100)

indices = find(sum(DF,2) ~= 0); % which time course differ
M = sqrt(mean([D(indices,:) - FilteredM(indices,:)].^2,2)); % data vs matlab filter
overfit = sum(M == 0) / size(M,1);
fprintf('matlab median filter matches exactly the data in %g%% of non match cases \n',overfit*100)

M = sqrt(mean([D(indices,:) - filteredM(indices,:)].^2,2)); % data vs spmup filter
goodfit = sum(M < 1) / size(M,1);
fprintf('spmup filter matches the data correctly in %g%% of cases (rms<1)\n',goodfit*100)

bad_fit = find(M>=1);
figure;
for i=1:length(bad_fit)
    subplot(1,2,1);
    plot([1:101],DataM(bad_fit(i),:),'*b','LineWidth',2);
    hold on; grid on; axis tight
    plot([1:101],DataM(bad_fit(i),:),'b','LineWidth',2);
    plot(filteredM(bad_fit(i),:),'r','LineWidth',2);
    title(['filtered data ' num2str(indices(i)) ' RMS ' num2str(M(bad_fit(i)))])
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
if exist('smooth.m','file')
    flags.method = 'moving';
    [despiked2,filtered2]  = spmup_despike(Data,[],flags);
    
    index = 1;
    for x=1:5
        for y=1:5
            for z=1:5
                filteredMM(index,:) = squeeze(filtered2(x,y,z,:));
                SS(index,:)  = squeeze(despiked(x,y,z,:));
                index = index+1;
            end
        end
    end
    
    DF = (filteredM(:,1:100) == filteredMM(:,1:100));
    DS = (S(:,1:100) == SS(:,1:100));
    
    figure
    subplot(1,2,1); plot(mean(filteredM),'LineWidth',3); hold on;
    plot(mean(filteredMM),'LineWidth',2); title('SPMUP filters median and moving avg'); axis tight; grid on
    subplot(1,2,2); plot(mean(S),'LineWidth',3); hold on;
    plot(mean(SS),'LineWidth',2); title('SPMUP data median and moving avg'); axis tight; grid on
end
