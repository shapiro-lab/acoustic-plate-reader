clear all

saveName = '/Users/Sanyo 1/Documents/MATLAB/Vantage-3.3.0-1710061400/Data/PlateReader/';
pathName = '/Users/Sanyo 1/Documents/MATLAB/Vantage-3.3.0-1710061400/Data/PlateReader/';
ExperimentDate = 'Raw';
SampleName = 'dC-Ser-S281+Ca';

pathName = fullfile(pathName, ExperimentDate, SampleName);
SubSet = GetSubDirs(pathName);
PlateCoordinate = strings(1,96);
total_n = length(SubSet);
PostCollapse = 1;
computeDiff = 1;
compar = [1 2];

for i = 1:total_n
    name_temp = char(SubSet(i));
    fileName = dir([pathName name_temp '/*.mat']);
    PlateCoordinate(i) = name_temp(end-2:end);
    Nf= length(fileName);
    j = 1;
    if i == 1
        h = waitbar(0,'Loading files');
        Im = cell(Nf,2,total_n);
        Imi = Im;
        voltage = nan(1,Nf);
        load([pathName char(SubSet(i)) '/' fileName(j).name]);
        voltage(j) = P.hv;
        Im{j,1,i} = ImgData.Imx;
        Im{j,2,i} = ImgData.Imb;
        X = cell(1,2);
        Z = cell(1,2);
        zmesh_temp = cell(1,2);
        xmesh_temp = cell(1,2);
        Z{1} = ImgData.Zx;
        X{1} = ImgData.Xx;
        Z{2} = ImgData.Zb;
        X{2} = ImgData.Xb;
        Xrange = [max(min(ImgData.Xx),min(ImgData.Xb)) min(max(ImgData.Xx),max(ImgData.Xb)) min(length(ImgData.Xx),length(ImgData.Xb))];
        Zrange = [max(min(ImgData.Zx),min(ImgData.Zb)) min(max(ImgData.Zx),max(ImgData.Zb)) min(length(ImgData.Zx),length(ImgData.Zb))];
        Xi = linspace(Xrange(1),Xrange(2),Xrange(3));
        Zi = linspace(Zrange(1),Zrange(2),Zrange(3));
        [Zi_mesh, Xi_mesh] = meshgrid(Xi,Zi);
        %Imi = nan(length(Zi),length(Xi),Nf, 2, total_n);
        if computeDiff
            %dImi = nan(length(Zi),length(Xi), 2, total_n);
            dImi = cell(Nf-1,2,total_n);
        end
        for k = 1:2
            [zmesh_temp{k},xmesh_temp{k}] = meshgrid(X{k},Z{k});
            Imi{j,k,i} = interp2(zmesh_temp{k},xmesh_temp{k}, Im{j,k,i},Zi_mesh,Xi_mesh);
        end
        waitbar((i*Nf+j)/(Nf*total_n))
    else
        load([pathName char(SubSet(i)) '/' fileName(j).name]);
        Im{j,1,i} = ImgData.Imx;
        Im{j,2,i} = ImgData.Imb;
        for k = 1:2
            Imi{j,k,i} = interp2(zmesh_temp{k},xmesh_temp{k}, Im{j,k,i},Zi_mesh,Xi_mesh);
        end
        waitbar((i*Nf+j)/(Nf*total_n))
    end
    for j = Nf:-1:2
        load([pathName char(SubSet(i)) '/' fileName(j).name]);
        Im{j,1,i} = ImgData.Imx;
        Im{j,2,i} = ImgData.Imb;
        voltage(j) = P.hv;
        for k = 1:2
            Imi{j,k,i} = interp2(zmesh_temp{k},xmesh_temp{k}, Im{j,k,i},Zi_mesh,Xi_mesh);
            if computeDiff && (j < Nf)
                dImi{j-1,k,i} = interp2(zmesh_temp{k},xmesh_temp{k}, (Im{j,k,i} - Im{Nf,k,i}),Zi_mesh,Xi_mesh);
            end
        end
        waitbar((i*Nf+(Nf+1-j))/(Nf*total_n))
    end
            
end
clear i j k;
close(h);
