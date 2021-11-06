% clear all
%close all

samplePath = fullfile(pathName, SampleName);
SubDirs = GetSubDirs(samplePath);
total_n = length(SubDirs); % total number of wells scanned

PlateCoordinate = strings(1,96);

% iterate through all wells
for well = 1:total_n
    name_temp = char(SubDirs(well));
    fileName = dir(fullfile(samplePath, name_temp, '*.mat'));
    PlateCoordinate(well) = name_temp(end-2:end); % get row and column of well
    Nf= length(fileName); % number of images per well
    
    % load in first image for all wells
    pressure = 1;
    % load in first image for first well
    if well == 1
        h = waitbar(0,'Loading your crappy data');
        Im = cell(Nf,2,total_n); % make Im with dimensions [# images per well, # imaging modes, # wells]
%         Imi = Im;
        
        
        load(fullfile(samplePath, char(SubDirs(well)), fileName(pressure).name));
        Im{pressure,1,well} = ImgData.Imx; % get xAM image
        Im{pressure,2,well} = ImgData.Imb; % get Bmode image
        
        % make cell arrays to store height and width vectors for xAM and
        % Bmode images
        X = cell(1,2);
        Z = cell(1,2);
        Z{1} = ImgData.Zx; % height vector for xAM image
        X{1} = ImgData.Xx; % width vector for xAM image
        Z{2} = ImgData.Zb; % height vector for Bmode image
        X{2} = ImgData.Xb; % width vector for Bmode image
        % in case xAM and Bmode images are not the same size, define X and
        % Z ranges such that they will occur in both images
        Xrange = [max(min(ImgData.Xx),min(ImgData.Xb)) min(max(ImgData.Xx),max(ImgData.Xb)) min(length(ImgData.Xx),length(ImgData.Xb))];
        Zrange = [max(min(ImgData.Zx),min(ImgData.Zb)) min(max(ImgData.Zx),max(ImgData.Zb)) min(length(ImgData.Zx),length(ImgData.Zb))];
        % make vectors for X and Z ranges
        Xi = linspace(Xrange(1),Xrange(2),Xrange(3));
        Zi = linspace(Zrange(1),Zrange(2),Zrange(3));
        % make meshgrid from X and Z ranges to use for interpolation
        xmesh_temp = cell(1,2);
        zmesh_temp = cell(1,2);
        [Zi_mesh, Xi_mesh] = meshgrid(Xi,Zi);
        % make 5D array to store interpolated images
        % Imi dimensions: zs, xs, pressures, imaging modes, wells
        Imi = nan(length(Zi),length(Xi),Nf, 2, total_n); 

        if computeDiff
            % make array to store pre-post-collapse difference images
            % dImi dimensions: zs, xs, imaging modes, wells
            dImi = nan(length(Zi),length(Xi), 2, total_n);
        end
        
        for mode = 1:2
            [zmesh_temp{mode},xmesh_temp{mode}] = meshgrid(X{mode},Z{mode});
            % linearly interpolate xAM and Bmode images with new
            % coordinates [Zi_mesh, Xi_mesh] to make all pixels the same
            % size
            Imi(:,:,pressure,mode,well) = interp2(zmesh_temp{mode},xmesh_temp{mode},Im{pressure,mode,well}, Zi_mesh,Xi_mesh);
        end
        
        waitbar((well*Nf+pressure)/(Nf*total_n))
    
    % load in first image for all other wells
    else
        load(fullfile(samplePath, char(SubDirs(well)), fileName(pressure).name));
        Im{pressure,1,well} = ImgData.Imx; % get xAM image
        Im{pressure,2,well} = ImgData.Imb; % get Bmode image
        
        for mode = 1:2
            % linearly interpolate xAM and Bmode images with new coordinates [Zi_mesh, Xi_mesh]
            Imi(:,:,pressure,mode,well) = interp2(zmesh_temp{mode},xmesh_temp{mode},Im{pressure,mode,well}, Zi_mesh,Xi_mesh);
        end
        
        waitbar((well*Nf+pressure)/(Nf*total_n))
        
    end
    
    % load in rest of images for all wells
    for pressure = 2:Nf
        load(fullfile(samplePath, char(SubDirs(well)), fileName(pressure).name));
        Im{pressure,1,well} = ImgData.Imx; % get xAM image
        Im{pressure,2,well} = ImgData.Imb; % get Bmode image
        
        for mode = 1:2
            % linearly interpolate xAM and Bmode images with new
            % coordinates [Zi_mesh, Xi_mesh] to make all pixels the same
            % size
            Imi(:,:,pressure,mode,well) = interp2(zmesh_temp{mode},xmesh_temp{mode},Im{pressure,mode,well}, Zi_mesh,Xi_mesh);
        end
        
        waitbar((well*Nf+pressure)/(Nf*total_n))
        
    end
    
    % get pre-post-collapse difference
    if computeDiff
        for mode = 1:2
            % linearly interpolate xAM and Bmode pre-post-collapse
            % difference images with new coordinates [Zi_mesh, Xi_mesh] to
            % make all pixels the same size
            dImi(:,:,mode,well) = interp2(zmesh_temp{mode},xmesh_temp{mode},(Im{compar(1),mode,well} - Im{compar(2),mode,well}), Zi_mesh,Xi_mesh);
        end
    end
            
end
close(h);

PlateSize = [P.xLines P.zLines]; % rows and columns of wells scanned
save(fullfile(samplePath, 'imgs.mat'), '-v7.3')

function [subDirsNames] = GetSubDirs(parentDir)
    % Get a list of all files and folders in folder.
    files    = dir(parentDir);
    names    = {files.name};
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir] & ...
              ~strcmp(names, '.') & ~strcmp(names, '..');
    % Extract only those that are directories.
    subDirsNames = string(names(dirFlags));
end
