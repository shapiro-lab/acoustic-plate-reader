clear all
close all

PlateSize = [8,12];
n_plates = 4;

% load data
data1 = load('/Volumes/GoogleDrive/.shortcut-targets-by-id/0B24ONICaZ0z9djczVE1ZR3BnWU0/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-lib-2/A-lib-T6A-A2/A-lib-T6A-A2_P1_R1-3_stable-37C_P_R2_C1_data_221103-13-24-37.mat');
data2 = load('/Volumes/GoogleDrive/.shortcut-targets-by-id/0B24ONICaZ0z9djczVE1ZR3BnWU0/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-lib-2/A-lib-T6A-A2/A-lib-T6A-A2_P2_R1-3_stable-37C_P_R1_C1_data_221110-11-44-32.mat');
data3 = load('/Volumes/GoogleDrive/.shortcut-targets-by-id/0B24ONICaZ0z9djczVE1ZR3BnWU0/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-lib-2/A-lib-T6A-A2/A-lib-T6A-A2_P3_R1-3_stable-37C_P_R1_C2_data_221102-12-25-13.mat');
data4 = load('/Volumes/GoogleDrive/.shortcut-targets-by-id/0B24ONICaZ0z9djczVE1ZR3BnWU0/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/A-lib-2/A-lib-T6A-A2/A-lib-T6A-A2_P4_R1-4_stable-37C_P_R1_C3_data_221110-12-23-00.mat');

%%
%combine arrays
CNRs = cat(4, squeeze(data1.sampCNR(:,:,:)), squeeze(data2.sampCNR(:,:,:)), squeeze(data3.sampCNR(:,:,:)), squeeze(data4.sampCNR(:,:,:)));
noiseROI_means = cat(4, squeeze(data1.noiseROI_means(:,:,:)), squeeze(data2.noiseROI_means(:,:,:)), squeeze(data3.noiseROI_means(:,:,:)), squeeze(data4.noiseROI_means(:,:,:)));
sampROI_means = cat(4, squeeze(data1.sampROI_means(:,:,:)), squeeze(data2.sampROI_means(:,:,:)), squeeze(data3.sampROI_means(:,:,:)), squeeze(data4.sampROI_means(:,:,:)));

max_CNRs = squeeze(max(CNRs, [], 1));
max_noiseROI_means = squeeze(max(noiseROI_means, [], 1));
max_sampROI_means = squeeze(max(sampROI_means, [], 1));

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