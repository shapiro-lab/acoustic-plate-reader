clear all
close all

% load data
data1 = load('/Volumes/GoogleDrive/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/B-lib-1/B-lib-1_plate4_rep1_stable-37C_P_1_4/B-lib-1_plate4_rep1_stable-37C_P_1_4_data_220404-16-18-03.mat');
data2 = load('/Volumes/GoogleDrive/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/B-lib-1/B-lib-1_plate4_rep2_stable-37C_P_1_2/220331_B-lib-1_plate4_rep2_stable-37C_P_1_2_data_220404-17-32-32.mat');
data3 = load('/Volumes/GoogleDrive/My Drive/Shapiro Lab Information/Data/Rob/96-well_plate_scans/GvpA-B-mutants/B-lib-1/B-lib-1_plate4_rep3_stable-37C_P_2_2/220331_B-lib-1_plate4_rep3_stable-37C_P_2_2_data_220404-17-46-39.mat');

%%
data = data1;

CNRs = cat(4, squeeze(data1.sampCNR(:,:,:)), squeeze(data2.sampCNR(:,:,:)), squeeze(data3.sampCNR(:,:,:)));
noiseROI_means = cat(4, squeeze(data1.noiseROI_mean(:,:,:)), squeeze(data2.noiseROI_mean(:,:,:)), squeeze(data3.noiseROI_mean(:,:,:)));
sampROIs = cat(4, squeeze(data1.sampROI(:,:,:)), squeeze(data2.sampROI(:,:,:)), squeeze(data3.sampROI(:,:,:)));

data.sampCNR = mean(CNRs,4);
data.noiseROI_mean = mean(noiseROI_means,4);
data.sampROI = mean(sampROIs,4);

data.sampCNR_stds = std(CNRs,0,4);
data.noiseROI_mean_stds = std(noiseROI_means,0,4);
data.sampROI_stds = std(sampROIs,0,4);

%%
save([data.saveName '_data_' datestr(now,'yymmdd-hh-MM-ss')],'-struct','data');






% %%
% split = 1; % number of voltage ramps
% NpS = 1; % number of replicates
% numcondition = 96; % number of unique samples
% PreV = 1:16; % Define pre-collapse voltage range
% % PostV = PreV(end)+1; % Define post-collapse voltage
% PostV = length(data1.voltage)-PreV(end)+1:length(data1.voltage); % Define post-collapse voltage range
% PreV_split = reshape(PreV,[],split)'; % split pre-collapse voltage range into individual ramps
% colors_group = linspecer(numcondition);
% %%
% figure; hold on;
% for samp = 1:numcondition
%     for ramp = 1:split
%         plot(data1.P.Vseq,data(PreV_split(ramp,:),samp),'LineWidth',3,'Color',[colors_group(samp,:) 1 - 0.3*(ramp-1)],'DisplayName',group_names(samp));
% %         pl = errorbar(vx,(mnormx_group(PreV_split(j,:),i)),(normxstd_group(PreV_split(j,:),i)),'LineWidth',3,'Color',[colors_group(i,:) 1 - 0.3*(j-1)],'DisplayName',group_names(i));
%     end
% end
% xlabel('Voltage (V)');
% ylabel(yLabel)
% legend;
% set(gca,'fontsize',16);