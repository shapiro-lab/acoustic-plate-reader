clear all
close all

PlateSize = [8,12];
n_plates = 4;

% load data
data1 = load('/Volumes/GoogleDrive-118305181097921812507/.shortcut-targets-by-id/0B24ONICaZ0z9djczVE1ZR3BnWU0/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-lib-1/A-lib_P1_R2-4_stable-37C_P_1_1_data_220727-11-54-47.mat');
data2 = load('/Volumes/GoogleDrive-118305181097921812507/.shortcut-targets-by-id/0B24ONICaZ0z9djczVE1ZR3BnWU0/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-lib-1/A-lib_P2_R1-3_stable-37C_data_220727-11-58-06.mat');
data3 = load('/Volumes/GoogleDrive-118305181097921812507/.shortcut-targets-by-id/0B24ONICaZ0z9djczVE1ZR3BnWU0/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-lib-1/A-lib_P3_R1-3_stable-37C_data_220727-12-02-59.mat');
data4 = load('/Volumes/GoogleDrive-118305181097921812507/.shortcut-targets-by-id/0B24ONICaZ0z9djczVE1ZR3BnWU0/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-lib-1/A-lib_P4_R1-3_stable-37C_P_2_4_data_220727-12-07-08.mat');

%% combine arrays
CNRs = cat(4, squeeze(data1.sampCNR(:,:,:)), squeeze(data2.sampCNR(:,:,:)), squeeze(data3.sampCNR(:,:,:)), squeeze(data4.sampCNR(:,:,:)));
noiseROI_means = cat(4, squeeze(data1.noiseROI_means(:,:,:)), squeeze(data2.noiseROI_means(:,:,:)), squeeze(data3.noiseROI_means(:,:,:)), squeeze(data4.noiseROI_means(:,:,:)));
noiseROI_stds = cat(4, squeeze(data1.noiseROI_stds(:,:,:)), squeeze(data2.noiseROI_stds(:,:,:)), squeeze(data3.noiseROI_stds(:,:,:)), squeeze(data4.noiseROI_stds(:,:,:)));
sampROI_means = cat(4, squeeze(data1.sampROI_means(:,:,:)), squeeze(data2.sampROI_means(:,:,:)), squeeze(data3.sampROI_means(:,:,:)), squeeze(data4.sampROI_means(:,:,:)));

%% calculate CNRs and SBRs
% calculate CNRs and SBRs for every well and pressure
sampCNR = (sampROI_means - noiseROI_means) ./ noiseROI_stds; % Calculate sample CNR
sampCNR_dB = 20 * log10(abs(sampCNR)); % Calculate sample CNR in dB
sampSBR = sampROI_means ./ noiseROI_means; % Calculate sample SBR
sampSBR_dB = 20 * log10(sampSBR); % Calculate sample SBR in dB

% calculate pre-/post-collapse difference and then calculate CNRs and SBRs
% from those values
sampCNR_diff = sampCNR(1:10,:,:,:) - sampCNR(11:20,:,:,:); % Calculate pre-/post-collapse difference CNR
sampCNR_diff_dB = 20 * log10(abs(sampCNR_diff)); % Calculate pre-/post-collapse difference CBR in dB
sampSBR_diff = sampSBR(1:10,:,:,:) - sampSBR(11:20,:,:,:); % Calculate pre-/post-collapse difference SBR
sampSBR_diff_dB = 20 * log10(abs(sampSBR_diff)); % Calculate pre-/post-collapse difference SBR in dB

%% take maxs and combine plates
max_CNRs = squeeze(max(CNRs, [], 1));

combined = table;
name1 = repmat([{'1'}, {'2'}, {'3'}, {'4'}],[96 1]);
name2 = repmat(data1.PlateCoordinate,[1 n_plates]);
combined.names = categorical(strcat(name1(:).', '_', name2)).';
combined.values = reshape(max_CNRs, [2,96*n_plates]).';
save([data1.saveName '_data-combined_' datestr(now,'yymmdd-hh-MM-ss') '.mat'],'combined');
combined_sorted = sortrows(combined,2,'descend');
writetable(combined_sorted,[data1.saveName '_quants.xlsx'])

%% plot samples from one plate sorted by CNR
figure;
y = squeeze(max(data1.sampCNR,[],1));
signals = table;
signals.names = categorical(data1.PlateCoordinate).';
signals.values = y(1,:).';
signals = sortrows(signals,2,'descend');
bar(reordercats(signals.names,cellstr(signals.names)), signals.values)
title('Max signal achieved at any voltage')
xlabel('Colony #');
ylabel('xAM CNR (dB)')

%% plot all samples
figure;
combined_sorted = sortrows(combined,2,'descend');
bar(reordercats(combined_sorted.names,cellstr(combined_sorted.names)), combined_sorted.values(:,1))
title('Max signal achieved at any voltage')
xlabel('Well');
ylabel('xAM CNR (dB)')
savefig([data1.saveName '_max-signal_combined.fig'])
pause(0.5)

figure;
bar(sort(combined.values(:,1),'descend'))
title('Max signal achieved at any voltage')
xlabel('Colony #');
ylabel('xAM CNR (dB)')
savefig([data1.saveName '_max-signal_combined.fig'])
pause(0.5)

figure;
histogram(combined.values(:,1),20)
title('Max signal achieved at any voltage')
xlabel('xAM CNR (dB)');
ylabel('# Colonies')
savefig([data1.saveName '_max-signal_combined_hist.fig'])