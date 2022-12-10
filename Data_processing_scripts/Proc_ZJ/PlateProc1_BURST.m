clear all

saveName = 'C:\Users\verasonics\Dropbox\GV Team\verasonics system\Vantage-4.6.2-RCH\Data\Wheesoo_221119_BURST/';
pathName = 'C:\Users\verasonics\Dropbox\GV Team\verasonics system\Vantage-4.6.2-RCH\Data\Wheesoo_221119_BURST/';
%ExperimentDate = '2022-08-06@19-51-20';
SampleName = 'GV3_2022-11-19@14-41-402';

pathName = fullfile(pathName, SampleName);
SubSet = GetSubDirs(pathName);
PlateCoordinate = strings(1,36);
total_n = length(SubSet);
disp_depth = [5 10];
BURST_filehead = append(extractBefore(SubSet(1),strlength(SubSet(1))-3),"_BURST");
for i = 1:total_n
    fileName = dir(fullfile(pathName,SubSet(i),BURST_filehead,'*.mat'));
    PlateCoordinate(i) = extractAfter(SubSet(i),strlength(SubSet(i))-3);
    Nf= length(fileName);
    for j = 1:Nf
        if i*j == 1
            h = waitbar(0,'Loading files');
            Im = cell(Nf,total_n);
            load(fullfile(pathName,SubSet(i),BURST_filehead,fileName(j).name));
            Xi = P.x;
            Zi = P.z;
            FrameCount = [1 P.numPreColFrames+1 P.numColFrames+P.numPreColFrames Nf];
            sum_temp = nan(length(FrameCount),total_n);
            count_temp = 1:length(FrameCount);
            iZd = find(Zi>=disp_depth(1)):find(Zi>=disp_depth(2));
            dImi = nan(length(Zi),length(Xi), 2, total_n);
            % iXd = find(Xi>=disp_width(1)):find(Xi>=disp_width(2));
            Im{j,i} = imData;
            sum_temp(j,i) = mean(mean(imData(iZd,:)));
        else
            load(fullfile(pathName,SubSet(i),BURST_filehead,fileName(j).name));
            Im{j,i} = imData;
            if ismember(j,FrameCount)
                sum_temp(count_temp(FrameCount==j),i) = mean(mean(imData(iZd,:)));
            end 
            waitbar((i*Nf+j)/(Nf*total_n))
        end
    end
    dImi(:,:,1,i) = Im{FrameCount(2),i} - Im{FrameCount(3),i};
    dImi(:,:,2,i) = Im{FrameCount(1),i} - Im{FrameCount(4),i};
end
clear i j;
close(h);
BURST_sum = reshape(sum_temp(2,:) - sum_temp(3,:),6,6)';

%BURST_sum = reshape(sum_temp(2,:) - sum_temp(3,:),12,8)
