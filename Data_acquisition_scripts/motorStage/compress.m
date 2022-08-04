clearvars;
clc;
disp('Finding Cavitation Files to Compress');
f = subdir('*Cavitation*.mat');

%%

fsizes = zeros(617,1);

s = zeros(617,1);

for i = 1:617
    s(i) = f(i).bytes;
end

[~,I] = sort(s);

f = f(I);
%%
clc
for i = 473:617
fn = f(i).name;

clear params
load(fn);

if exist('params');
a = whos('params');
disp(sprintf('(%03.0f) %s (%1.2f MB)', i,fn(find(fn == '\',1,'last')+9:end-4), a.bytes / 1e6));

flag = 0;

if isfield(params, 'ScriptCode')
if isfield(params, 'Results')
    if isfield(params, 'Scope')
        disp('- Removing Scope File');
        params = rmfield(params, 'Scope');
        flag = 1;
    end
    
    if isfield(params.Results, 'w')
        disp('- Removing w');
        params.Results = rmfield(params.Results, 'w');
        flag = 1;
    end
       
    if isfield(params.Results, 't_sample')
        disp('- Removing t_sample');
        params.Results = rmfield(params.Results, 't_sample');
        flag = 1;
    end
    
    if isfield(params.Results, 'FFT')
        disp('- Removing FFT');
        params.Results = rmfield(params.Results, 'FFT');
        flag = 1;
    end
    
    if isfield(params.Results, 'waveforms')
        if isfield(params.Results, 'data')
            params.Results = rmfield(params.Results, 'data');
            disp('- Removing data')
            flag = 1;
        end
    end
    
    if isfield(params.Results, 'waveform')
    if isfield(params.Results, 'data')
        params.Results = rmfield(params.Results, 'data');
        disp('- Removing data')
        flag = 1;
    end
    if ~isfield(params.Results.waveform, 'XDataComp')
        disp('- Compressing waveform');
               
        pre = params.Results.waveform;        
        params.Results.waveform = sub_Data_CompressWaveform(params.Results.waveform);
        post = sub_Data_DecompressWaveform(params.Results.waveform);
        
        if     ~isequal(pre.XData, post.XData) && ~isequal(pre.XData, post.XData');
            subplot(311); plot(pre.XData, pre.YData); legend('Pre-Compression')
            subplot(312); plot(post.XData, post.YData); legend('Post-Compression')
            if size(post.XData) == size(pre.XData);
            subplot(313); plot(post.XData - pre.XData); legend('Post X - Pre X')
                if max(abs(post.XData - pre.XData)) > post.XData(2) - post.XData(1)
                    disp('- Significant X Data error, press any key to continue')
                    pause
                else
                    disp('- Insignificant X Data error, compression fine')
                    flag = 1;
                end
            else
            subplot(313); plot(post.XData' - pre.XData); legend('Post X - Pre X')
                if max(abs(post.XData' - pre.XData)) > post.XData(2) - post.XData(1)
                    disp('- Significant X Data error, press any key to continue')
                    pause
                else
                    disp('- Insignificant X Data error, compression fine')
                    flag = 1;
                end
            end

        elseif ~isequal(pre.YData, post.YData) && ~isequal(pre.YData, post.YData');
            subplot(311); plot(pre.XData, pre.YData); legend('Pre-Compression')
            subplot(312); plot(post.XData, post.YData); legend('Post-Compression')
            if size(post.XData) == size(pre.XData);
            subplot(313); plot(post.XData - pre.XData); legend('Post X - Pre X')
            else
            subplot(313); plot(post.XData' - pre.XData); legend('Post X - Pre X')
            end
            disp('- Y Data error, how did that happen!  Press any key to continue')
            pause
        else
            if ~isfield(params.Results.waveform, 'DidNotCompress')
            flag = 1;
            end
        end
        
        
    end
    end
    
    if flag
       a = whos('params');
       disp(sprintf('- Compressed Size (%1.2f MB)', a.bytes / 1e6));
       save(fn, 'params', '-v7.3');
       clf;
    end
end
end
end

end