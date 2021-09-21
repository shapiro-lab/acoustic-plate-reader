close all
fft_pts = numel(t);
w = (0:fft_pts-1)./fft_pts.*fs;

hpFilt = designfilt('highpassiir', 'FilterOrder', 8, ...
         'PassbandFrequency', 1e6, 'PassbandRipple', 0.2,...
         'SampleRate', fs);

B = filtfilt(hpFilt,A);
     
subplot(311); plot(t, A, 'r', t, B, 'b')
subplot(312); plot(w, abs(fft(A)), 'r'); xlim([0 max(w)/2]); hold on;

for i = 1:20
    plot(494e3 * i, 0, 'k^')
end

subplot(313); plot(w, abs(fft(B)), 'b'); xlim([0 max(w)/2]); hold on;
for i = 1:20
    plot(494e3 * i, 0, 'k^')
end
