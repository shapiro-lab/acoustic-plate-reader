   %overlay = 1;
showf5 = 0;
xROI_size = 10;
zROI_size = 40;
zoffset = 80;
testfilter = 1;
rangeoffset = 10;
for testInd = 1:40
    %Ind = 59;
    Imt = Im{1,2,testInd};
    
    %dImi(:,:,imgMode,Ind);
    filter_size1 = [20 5];
    filter_size2 = [5 5];
    filter_sigma = 2.5;
    filter_box = [9 1];
    %SizeThreshold = 300;

    threshold_ratio = 1;
    noise_stdratio = 1;
    auto_noise = 1;
    noiseSize = [2 2];
    sample_depth = [3 8];
    noise_slice = [2 3];
    sample_width = 10;
    iZt = find(Zi>=sample_depth(1)):find(Zi>=sample_depth(2));
    iZd_noise = find(Zi>=noise_slice(1)):find(Zi>=noise_slice(2));
    filter_depth = zeros(size(Imt));
    filter_depth(iZt,:) = 1;
    Zit = Zi(iZt);

    if auto_noise
        noiseROIt = mean(mean(Imt(iZd_noise,:))) + noise_stdratio * std2(Imt(iZd_noise,:));
    else
        noiseROIt = noiseROI(:,Ind);
    end

    Imt_f1 = medfilt2(Imt(iZt,:),filter_size1);
    threshold1 = noiseROIt*1;
    Imt_f2 = (Imt_f1 > threshold1);
    if testfilter
        figure;
        imagesc(Xi, Zit, Imt_f2); axis image;
    end
    
    Xm = repmat(1:length(Xi),length(Zit),1).* Imt_f2;
    %Zm = repmat((1:length(Zit))',1,length(Xi)).* Imt_f2;
    bx1 = max(Xm,[],2) - min(Xm,[],2);
    bx2 = sum(Imt_f2,2);
    [~,ztop] = max(bx1);
    izbot = max(any(Imt_f2,2).*(1:length(Zit))');
    izbot_range = -ztop + izbot - rangeoffset;
    ixcenter = round(sum(sum((Xm(izbot - izbot_range:izbot,:)))) / nnz(Xm(izbot - izbot_range:izbot,:)));
    
    
    
    
    
%     b = bwboundaries(Imt_f2,4);
%     b_size = cell2mat(cellfun(@size,b,'UniformOutput',false));
%     [~,bn] = max(b_size(:,1));
%     boundary = b{bn};
% 
%     bx = unique(boundary(:,2));
%     bz = unique(boundary(:,1));
%     bz1 = nan(3,length(bx));
%     bx1 = nan(3,length(bz));
%     for i = 1:length(bx)
%         temp = boundary(boundary(:,2)==bx(i),1);
%         bz1(1,i) = bx(i);
%         bz1(2,i) = length(temp);
%         bz1(3,i) = max(temp) - min(temp);
%     end
% 
%     [~,I_temp] = sort(bz1(3,:),'descend');
%     bz1 = bz1(:,I_temp);
%     
%     xtemp = [];
%     for i = 1:length(bz)
%         temp = boundary(boundary(:,1)==bz(i),2);
%         bx1(1,i) = bz(i);
%         bx1(2,i) = length(temp);
%         bx1(3,i) = max(temp) - min(temp);
%         if i > length(bz)-10
%             xtemp = [xtemp;unique(temp)];
%         end
%     end
%     xtemp = round(median(xtemp));
%     [~,I_temp] = sort(bx1(3,:),'descend');
%     bx1 = bx1(:,I_temp);
% 
    Ind_dX = [ixcenter-xROI_size ixcenter+xROI_size];
% 
%     %temp = bx1(1,bx1(3,:)==(Ind_dX0(2)-Ind_dX0(1)));
    Ind_dZ = [ztop-zROI_size ztop+zROI_size];
    Ind_dZ = Ind_dZ + zoffset;
    %Ind_dX = [Ind_dX0(1) + xoffset Ind_dX0(2) - xoffset];

    figure;
    imagesc(Xi, Zit, 20*log10(abs(Imt(iZt,:))), [20 80]); axis image;colormap hot
    hold on;
    rectangle('Position',[Xi(Ind_dX(1)) Zit(Ind_dZ(1)) Xi(Ind_dX(2))-Xi(Ind_dX(1)) Zit(Ind_dZ(2))-Zit(Ind_dZ(1))],'EdgeColor','w','LineWidth',2)
    pause(1);
    %close all;
end
% Imt2 = zeros(size(Imt_f1));
% for i = 1:length(boundary(:,1))
% Imt2(boundary(i,1),boundary(i,2)) = 1;
% end
% 
% figure;
% imagesc(Xi, Zit, Imt2, [0 1]);axis image;colormap hot