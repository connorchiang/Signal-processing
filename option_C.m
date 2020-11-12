load('segAllData.mat');

% DFT for 4 mats, Number of samples in signal frame N = 1024
DFT_Orig_phAA = abs(fft(segOrig_phAA,1024,2));
DFT_Orig_phS = abs(fft(segOrig_phS,1024,2));
DFT_Filt_phAA = abs(fft(segFilt_phAA,1024,2));
DFT_Filt_phS = abs(fft(segFilt_phS,1024,2));

% DFT converts to (dB)
DFT_Orig_phAA_dB = 20*log10(DFT_Orig_phAA);
DFT_Orig_phS_dB = 20*log10(DFT_Orig_phS);
DFT_Filt_phAA_dB = 20*log10(DFT_Filt_phAA);
DFT_Filt_phS_dB = 20*log10(DFT_Filt_phS);

figure(1);
subplot(2,2,1);
% Fs=8000Hz, N=1024
% plot magnitude spectrum for 
% the first occurrence of the phoneme ¡®s¡¯ and¡®aa' in the files:¡®wavOrig/MDPK0/SA1.wav¡¯ and ¡®wavFilt/MDPK0/SA1.wav¡¯.
plot([0:511]*8000/1024,DFT_Orig_phAA_dB(24,1:512)); 
grid on; xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); 
title('magnitude spectrum (first aa in wavOrig/MDPK0/SA1.wav)'); axis([0 4000 -60 40]);
subplot(2,2,2);
plot([0:511]*8000/1024,DFT_Orig_phS_dB(28,1:512));
grid on; xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); 
title('magnitude spectrum (first s in wavOrig/MDPK0/SA1.wav)'); axis([0 4000 -60 40]);
subplot(2,2,3);
plot([0:511]*8000/1024,DFT_Filt_phAA_dB(24,1:512));
grid on; xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); 
title('magnitude spectrum (first aa in wavFilt/MDPK0/SA1.wav)'); axis([0 4000 -60 40]);
subplot(2,2,4);
plot([0:511]*8000/1024,DFT_Filt_phS_dB(28,1:512));
grid on; xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); 
title('magnitude spectrum (first s in wavFilt/MDPK0/SA1.wav)'); axis([0 4000 -60 40]);

% 4 arrays which store average energy 
% first column stores average energy for low frequency (0-2k Hz)
% second column for high frequency (2k-4k Hz)

enLF_andHF_orig_phAA_0 = zeros(size(DFT_Orig_phAA,1),2);
enLF_andHF_orig_phS_0 = zeros(size(DFT_Orig_phS,1),2);
enLF_andHF_filt_phAA_0 = zeros(size(DFT_Filt_phAA,1),2);
enLF_andHF_filt_phS_0 = zeros(size(DFT_Filt_phS,1),2);

% Caculate average energy for phoneme AA
for i=1:1:size(DFT_Orig_phAA,1)
    % Low frequency (0-2k HZ)
    for j=1:1:256
        enLF_andHF_orig_phAA_0(i,1) = enLF_andHF_orig_phAA_0(i,1) + DFT_Orig_phAA(i,j)^2;
        enLF_andHF_filt_phAA_0(i,1) = enLF_andHF_filt_phAA_0(i,1) + DFT_Filt_phAA(i,j)^2;
    end
    
    %High frequency (2k-4k Hz)
    for j=257:1:512
        enLF_andHF_orig_phAA_0(i,2) = enLF_andHF_orig_phAA_0(i,2) + DFT_Orig_phAA(i,j)^2;
        enLF_andHF_filt_phAA_0(i,2) = enLF_andHF_filt_phAA_0(i,2) + DFT_Filt_phAA(i,j)^2;
    end
end
% aveEn = average{|X(k)|.^2}_k
enLF_andHF_orig_phAA_0 = enLF_andHF_orig_phAA_0 / 256;
enLF_andHF_filt_phAA_0 = enLF_andHF_filt_phAA_0 / 256;

% Caculate average energy for phoneme S
for i=1:1:size(DFT_Orig_phS,1)
    % Low frequency (0-2k HZ)
    for j=1:1:256
        enLF_andHF_orig_phS_0(i,1) = enLF_andHF_orig_phS_0(i,1) + DFT_Orig_phS(i,j)^2;
        enLF_andHF_filt_phS_0(i,1) = enLF_andHF_filt_phS_0(i,1) + DFT_Filt_phS(i,j)^2;
    end
    %High frequency (2k-4k Hz)
    for j=257:1:512
        enLF_andHF_orig_phS_0(i,2) = enLF_andHF_orig_phS_0(i,2) + DFT_Orig_phS(i,j)^2;
        enLF_andHF_filt_phS_0(i,2) = enLF_andHF_filt_phS_0(i,2) + DFT_Filt_phS(i,j)^2;
    end
end
% aveEn = average{|X(k)|.^2}_k
enLF_andHF_orig_phS_0 = enLF_andHF_orig_phS_0 / 256;
enLF_andHF_filt_phS_0 = enLF_andHF_filt_phS_0 / 256;

% Energy converts to (dB)
enLF_andHF_orig_phAA = 20*log10(enLF_andHF_orig_phAA_0);
enLF_andHF_filt_phAA = 20*log10(enLF_andHF_filt_phAA_0);
enLF_andHF_orig_phS = 20*log10(enLF_andHF_orig_phS_0);
enLF_andHF_filt_phS = 20*log10(enLF_andHF_filt_phS_0);

figure(2);
subplot(2,4,1);
histogram(enLF_andHF_orig_phAA(:,1),-80:10:40); ylabel('number'), xlabel('average energy (dB)') ; 
title('en LF orig phAA'); axis([-80 40 0 40]); grid on;
subplot(2,4,2);
histogram(enLF_andHF_orig_phAA(:,2),-80:10:40); ylabel('number'), xlabel('average energy (dB)') ; 
title('en HF orig phAA'); axis([-80 40 0 40]); grid on;
subplot(2,4,3);
histogram(enLF_andHF_orig_phS(:,1),-80:10:40); ylabel('number'), xlabel('average energy (dB)') ; 
title('en LF orig phS'); axis([-80 40 0 40]); grid on;
subplot(2,4,4);
histogram(enLF_andHF_orig_phS(:,2),-80:10:40); ylabel('number'), xlabel('average energy (dB)') ; 
title('en HF orig phS'); axis([-80 40 0 40]); grid on;
subplot(2,4,5);
histogram(enLF_andHF_filt_phAA(:,1),-80:10:40); ylabel('number'), xlabel('average energy (dB)') ; 
title('en LF filt phAA'); axis([-80 40 0 40]); grid on;
subplot(2,4,6);
histogram(enLF_andHF_filt_phAA(:,2),-80:10:40); ylabel('number'), xlabel('average energy (dB)') ; 
title('en HF filt phAA'); axis([-80 40 0 40]); grid on;
subplot(2,4,7);
histogram(enLF_andHF_filt_phS(:,1),-80:10:60); ylabel('number'), xlabel('average energy (dB)') ; 
title('en LF filt phS'); axis([-80 40 0 40]); grid on;
subplot(2,4,8);
histogram(enLF_andHF_filt_phS(:,2),-80:10:60); ylabel('number'), xlabel('average energy (dB)') ; 
title('en HF filt phS'); axis([-80 40 0 40]); grid on;