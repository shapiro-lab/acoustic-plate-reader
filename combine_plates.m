clear all
close all

PlateSize = [8,12];

% load data
data1 = load('/Volumes/GoogleDrive/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/B-lib-1_old/B-lib-1_plate1_rep1-3_stable-37C_data_220406-16-41-01.mat');
data2 = load('/Volumes/GoogleDrive/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/B-lib-1_old/B-lib-1_plate2_rep1-3_stable-37C_data_220406-16-46-56.mat');
data3 = load('/Volumes/GoogleDrive/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/B-lib-1_old/B-lib-1_plate3_rep1-3_stable-37C_data_220406-16-49-06.mat');
data4 = load('/Volumes/GoogleDrive/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/B-lib-1_old/B-lib-1_plate4_rep1-3_stable-37C_data_220406-16-51-04.mat');

%%
%combine arrays
CNRs = cat(4, squeeze(data1.sampCNR(:,:,:)), squeeze(data2.sampCNR(:,:,:)), squeeze(data3.sampCNR(:,:,:)), squeeze(data4.sampCNR(:,:,:)));
noiseROI_means = cat(4, squeeze(data1.noiseROI_means(:,:,:)), squeeze(data2.noiseROI_means(:,:,:)), squeeze(data3.noiseROI_means(:,:,:)), squeeze(data4.noiseROI_means(:,:,:)));
sampROI_means = cat(4, squeeze(data1.sampROI_means(:,:,:)), squeeze(data2.sampROI_means(:,:,:)), squeeze(data3.sampROI_means(:,:,:)), squeeze(data4.sampROI_means(:,:,:)));

max_CNRs = squeeze(max(CNRs, [], 1));
max_noiseROI_means = squeeze(max(noiseROI_means, [], 1));
max_sampROI_means = squeeze(max(sampROI_means, [], 1));

combined = reshape(max_CNRs, [2,384]);

bar(sort(combined(1,:),'descend'))
title('Max signal achieved at any voltage')
xlabel('Colony');
ylabel('xAM CNR (dB)')