clear all
close all

PlateSize = [8,12];

% load data
data1 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-lib-1/A-lib_P1_R2-4_stable-37C_P_1_1_data_220727-11-54-47.mat');
data2 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-lib-1/A-lib_P2_R1-3_stable-37C_data_220727-11-58-06.mat');
data3 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-lib-1/A-lib_P3_R1-3_stable-37C_data_220727-12-02-59.mat');
data4 = load('/Users/Rob/Library/CloudStorage/GoogleDrive-rchurt@caltech.edu/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-lib-1/A-lib_P4_R1-3_stable-37C_P_2_4_data_220727-12-07-08.mat');

%%
%combine arrays
CNRs = cat(4, squeeze(data1.sampCNR(:,:,:)), squeeze(data2.sampCNR(:,:,:)), squeeze(data3.sampCNR(:,:,:)), squeeze(data4.sampCNR(:,:,:)));
noiseROI_means = cat(4, squeeze(data1.noiseROI_means(:,:,:)), squeeze(data2.noiseROI_means(:,:,:)), squeeze(data3.noiseROI_means(:,:,:)), squeeze(data4.noiseROI_means(:,:,:)));
sampROI_means = cat(4, squeeze(data1.sampROI_means(:,:,:)), squeeze(data2.sampROI_means(:,:,:)), squeeze(data3.sampROI_means(:,:,:)), squeeze(data4.sampROI_means(:,:,:)));

max_CNRs = squeeze(max(CNRs, [], 1));
max_noiseROI_means = squeeze(max(noiseROI_means, [], 1));
max_sampROI_means = squeeze(max(sampROI_means, [], 1));

combined = reshape(max_CNRs, [2,384]);
save([data1.saveName '_data-combined_' datestr(now,'yymmdd-hh-MM-ss')],'combined');


%% plot samples from one plate sorted by CNR
figure;
y = squeeze(max(data1.sampCNR,[],1));
s=table;
s.names = categorical(data1.PlateCoordinate).';
s.values = y(1,:).';
s = sortrows(s,2,'descend');
bar(reordercats(s.names,cellstr(s.names)), s.values)
title('Max signal achieved at any voltage')
xlabel('Colony');
ylabel('xAM CNR (dB)')

%% plot all samples
figure;
bar(sort(combined(1,:),'descend'))
title('Max signal achieved at any voltage')
xlabel('Colony');
ylabel('xAM CNR (dB)')
savefig([data1.saveName '_max-signal_combined'])
pause(0.5)

figure;
histogram(combined(1,:),20)
title('Max signal achieved at any voltage')
xlabel('xAM CNR (dB)');
ylabel('# Colonies')
savefig([data1.saveName '_max-signal_combined_hist'])