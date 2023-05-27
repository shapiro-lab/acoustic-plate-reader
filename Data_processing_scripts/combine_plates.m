clear all
close all

PlateSize = [8,12];
n_plates = 4;

% load data
<<<<<<< Updated upstream
data1 = load('G:\.shortcut-targets-by-id\0B24ONICaZ0z9djczVE1ZR3BnWU0\Shapiro Lab Information\Data\Rob\APR_scans\GvpA-B-mutants\A-lib-2\A-lib-K22R-A2\A-lib-K22R-A2_P1_R1-3_stable-37C_quants_230505-19-21-04.mat');
data2 = load('G:\.shortcut-targets-by-id\0B24ONICaZ0z9djczVE1ZR3BnWU0\Shapiro Lab Information\Data\Rob\APR_scans\GvpA-B-mutants\A-lib-2\A-lib-K22R-A2\A-lib-K22R-A2_P2_R1-3_stable-37C_quants_230505-19-22-20.mat');
data3 = load('G:\.shortcut-targets-by-id\0B24ONICaZ0z9djczVE1ZR3BnWU0\Shapiro Lab Information\Data\Rob\APR_scans\GvpA-B-mutants\A-lib-2\A-lib-K22R-A2\A-lib-K22R-A2_P3_R1-3_stable-37C_quants_230505-19-23-40.mat');
data4 = load('G:\.shortcut-targets-by-id\0B24ONICaZ0z9djczVE1ZR3BnWU0\Shapiro Lab Information\Data\Rob\APR_scans\GvpA-B-mutants\A-lib-2\A-lib-K22R-A2\A-lib-K22R-A2_P4_R1-3_stable-37C_quants_230505-19-24-47.mat');

%% combine arrays
%make structure to hold data
quants = data2;
=======
data1 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/B-lib-1/B-lib_P1_R134_stable-37C_P_1_1_data_220727-11-03-28.mat');
data2 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/B-lib-1/B-lib_P2_R134_stable-37C_P_1_2_data_220727-10-57-12.mat');
data3 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/B-lib-1/B-lib_P3_R134_stable-37C_P_1_3_data_220727-10-44-56.mat');
data4 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/B-lib-1/B-lib_P4_R2-4_stable-37C_P_1_4_data_220727-10-51-17.mat');
>>>>>>> Stashed changes

%combine arrays
quants.CNR_diffs = cat(4, squeeze(data1.sampCNR_diff(:,:,:)), squeeze(data2.sampCNR_diff(:,:,:)), squeeze(data3.sampCNR_diff(:,:,:)), squeeze(data4.sampCNR_diff(:,:,:)));
quants.SBR_diffs = cat(4, squeeze(data1.sampSBR_diff(:,:,:)), squeeze(data2.sampSBR_diff(:,:,:)), squeeze(data3.sampSBR_diff(:,:,:)), squeeze(data4.sampSBR_diff(:,:,:)));
quants.noiseROI_means = cat(4, squeeze(data1.noiseROI_means(:,:,:)), squeeze(data2.noiseROI_means(:,:,:)), squeeze(data3.noiseROI_means(:,:,:)), squeeze(data4.noiseROI_means(:,:,:)));
quants.sampROI_means = cat(4, squeeze(data1.sampROI_means(:,:,:)), squeeze(data2.sampROI_means(:,:,:)), squeeze(data3.sampROI_means(:,:,:)), squeeze(data4.sampROI_means(:,:,:)));

%% take maxs and combine plates
max_SBRs = squeeze(max(quants.SBR_diffs, [], 1));

combined = table;
name1 = repmat([{'1'}, {'2'}, {'3'}, {'4'}],[96 1]);
name2 = repmat(data1.PlateCoordinate,[1 n_plates]);
combined.names = categorical(strcat(name1(:).', '_', name2)).';
combined.values = reshape(max_SBRs, [2,96*n_plates]).';
save([data1.saveName '_data-combined_' datestr(now,'yymmdd-hh-MM-ss') '.mat'],'combined');
combined_sorted = sortrows(combined,2,'descend');
writetable(combined_sorted,[data1.saveName '_quants.xlsx'])

%% plot samples from one plate sorted by CNR
figure;
y = squeeze(max(data1.sampSBR_diff,[],1));
signals = table;
signals.names = categorical(data1.PlateCoordinate).';
signals.values = y(1,:).';
signals = sortrows(signals,2,'descend');
bar(reordercats(signals.names,cellstr(signals.names)), signals.values)
title('Max signal achieved at any voltage')
xlabel('Colony #');
ylabel('xAM difference SBR')

%% plot all samples
figure;
combined_sorted = sortrows(combined,2,'descend');
bar(reordercats(combined_sorted.names,cellstr(combined_sorted.names)), combined_sorted.values(:,1))
title('Max signal achieved at any voltage')
xlabel('Well');
ylabel('xAM difference SBR')
savefig([data1.saveName '_max-signal_combined.fig'])
pause(0.5)

figure;
bar(sort(combined.values(:,1),'descend'))
title('Max signal achieved at any voltage')
xlabel('Colony #');
ylabel('xAM difference SBR')
savefig([data1.saveName '_max-signal_combined.fig'])
pause(0.5)

figure;
histogram(combined.values(:,1),20)
title('Max signal achieved at any voltage')
xlabel('xAM difference SBR');
ylabel('# Colonies')
savefig([data1.saveName '_max-signal_combined_hist.fig'])