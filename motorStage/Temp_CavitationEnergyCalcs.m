clc
clearvars

files = subdir('*Cavitation*.mat');
disp(sprintf('Total files found: %1.0f', numel(files)));
disp(' ');
pause(3);

%%

energycalcs = [];

for i = 1:numel(files);
    disp(sprintf('(%02.1f%%) %s', i / numel(files)*100, files(i).name))
        

    clear params
    disp(' - Loading')
    success = 1;
    try; 
        load(files(i).name, 'params'); 
    catch; 
        disp(' - Could not load params')
        success = 0;
    end;
    
    if success
    if exist('params', 'var')
    if isfield(params, 'Results')
        if isfield(params, 'ScriptCode')


            eng = 0;
            n_wvs = 0;
            
            if isfield(params.Results, 'waveforms')
                n_wvs = numel(params.Results.waveforms);
                
                if iscell(params.Results.waveforms)
                for j = 1:numel(params.Results.waveforms)
                    y = params.Results.waveforms{j}.YData(:);
                    y = y - mean(y);
                    eng = eng +  sum(y.^ 2);
                end
                else
                for j = 1:numel(params.Results.waveforms)
                    y = params.Results.waveforms(j).YData(:);
                    y = y - mean(y);
                    eng = eng +  sum(y.^ 2);
                end
                end
                
            elseif isfield(params.Results, 'waveform')
                n_wvs = numel(params.Results.waveform);
                
                if iscell(params.Results.waveform)
                for j = 1:numel(params.Results.waveform)
                    y = params.Results.waveform{j}.YData(:);
                    y = y - mean(y);
                    eng = eng +  sum(y.^ 2);
                end
                else
                for j = 1:numel(params.Results.waveform)
                    y = params.Results.waveform(j).YData(:);
                    y = y - mean(y);
                    eng = eng +  sum(y.^ 2);
                end
                end
            
            else
                error('boop')
            end
            
            disp(sprintf(' - %1.3e', eng))
            
            ei = numel(energycalcs) + 1;
            energycalcs{ei}.name = files(i).name;
            energycalcs{ei}.energy = eng;
            energycalcs{ei}.params = rmfield(params, 'Results');
            energycalcs{ei}.n = n_wvs;
            
        end
    end
    end
    end
end

%%
clear params;
save('cavitation_energies.mat', 'energycalcs');

%%
clc
for i = 1:numel(energycalcs)
    str = energycalcs{i}.name;
    str = str(find(str == '\', 1, 'last')+9:end-4);
    disp(sprintf('%s, %1.3e, %1.0f, %1.0f', str, energycalcs{i}.energy, energycalcs{i}.params.SG.Waveform.frequency, energycalcs{i}.n))
end


%%
fileID = fopen('cavitation_energies.dat', 'w');
for i = 1:numel(energycalcs)
    str = energycalcs{i}.name;
    str = str(find(str == '\', 1, 'last')+9:end-4);
    fprintf(fileID,'%s, %1.3e, %1.0f, %1.0f\n', str, energycalcs{i}.energy, energycalcs{i}.params.SG.Waveform.frequency,  energycalcs{i}.n);
end
fclose(fileID);
