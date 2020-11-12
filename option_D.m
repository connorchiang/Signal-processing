load('enLF_andHF_orig_phAA.mat');
load('enLF_andHF_orig_phS.mat');
load('enLF_andHF_filt_phAA.mat');
load('enLF_andHF_filt_phS.mat');

% Expectation
mu_LF_andHF_orig_phAA = mean(enLF_andHF_orig_phAA);
mu_LF_andHF_filt_phAA = mean(enLF_andHF_filt_phAA);
mu_LF_andHF_orig_phS = mean(enLF_andHF_orig_phS);
mu_LF_andHF_filt_phS = mean(enLF_andHF_filt_phS);

mu_all = [mu_LF_andHF_orig_phAA; mu_LF_andHF_filt_phAA; mu_LF_andHF_orig_phS; mu_LF_andHF_filt_phS];

% Varience
var_LF_andHF_orig_phAA = var(enLF_andHF_orig_phAA,1,1);
var_LF_andHF_filt_phAA = var(enLF_andHF_filt_phAA,1,1);
var_LF_andHF_orig_phS = var(enLF_andHF_orig_phS,1,1);
var_LF_andHF_filt_phS = var(enLF_andHF_filt_phS,1,1);

var_all = [var_LF_andHF_orig_phAA; var_LF_andHF_filt_phAA; var_LF_andHF_orig_phS; var_LF_andHF_filt_phS];

% A table with the values of the parameters of the Gaussian PDFs modelling each data
Gaussian_Table = table (mu_all(:,1), mu_all(:,2), var_all(:,1),var_all(:,2), ...
            'VariableNames',{'mu_LF','mu_HF','var_LF','var_HF'},...
            'RowNames',{'orig_phAA';'filt_phAA';'orig_phS';'filt_phS'});

% plot Gaussian PDF
p_LF_orig_phAA = gaussPDF(enLF_andHF_orig_phAA(:,1), mu_LF_andHF_orig_phAA(1),var_LF_andHF_orig_phAA(1));
p_HF_orig_phAA = gaussPDF(enLF_andHF_orig_phAA(:,2), mu_LF_andHF_orig_phAA(2),var_LF_andHF_orig_phAA(2));
p_LF_orig_phS = gaussPDF(enLF_andHF_orig_phS(:,1), mu_LF_andHF_orig_phS(1),var_LF_andHF_orig_phS(1));
p_HF_orig_phS = gaussPDF(enLF_andHF_orig_phS(:,2), mu_LF_andHF_orig_phS(2),var_LF_andHF_orig_phS(2));
p_LF_filt_phAA = gaussPDF(enLF_andHF_filt_phAA(:,1), mu_LF_andHF_filt_phAA(1),var_LF_andHF_filt_phAA(1));
p_HF_filt_phAA = gaussPDF(enLF_andHF_filt_phAA(:,2), mu_LF_andHF_filt_phAA(2),var_LF_andHF_filt_phAA(2));
p_LF_filt_phS = gaussPDF(enLF_andHF_filt_phS(:,1), mu_LF_andHF_filt_phS(1),var_LF_andHF_filt_phS(1));
p_HF_filt_phS = gaussPDF(enLF_andHF_filt_phS(:,2), mu_LF_andHF_filt_phS(2),var_LF_andHF_filt_phS(2));

subplot(2,4,1);
plot(enLF_andHF_orig_phAA(:,1),p_LF_orig_phAA, '.');
xlabel('average energy (dB)'), ylabel('probability') ; 
title('GaussPDF for (en LF orig phAA)'); axis([-80 40 0 0.05]); grid on;
subplot(2,4,2);
plot(enLF_andHF_orig_phAA(:,2),p_HF_orig_phAA, '.');
xlabel('average energy (dB)'), ylabel('probability') ; 
title('GaussPDF for (en HF orig phAA)');  axis([-80 40 0 0.05]); grid on;
subplot(2,4,3);
plot(enLF_andHF_orig_phS(:,1),p_LF_orig_phS, '.');
xlabel('average energy (dB)'), ylabel('probability') ; 
title('GaussPDF for (en LF orig phS)'); axis([-80 40 0 0.05]); grid on;
subplot(2,4,4);
plot(enLF_andHF_orig_phS(:,2),p_HF_orig_phS, '.');
xlabel('average energy (dB)'), ylabel('probability') ; 
title('GaussPDF for (en HF orig phS)'); axis([-80 40 0 0.05]); grid on;
subplot(2,4,5);
plot(enLF_andHF_filt_phAA(:,1),p_LF_filt_phAA, '.');
xlabel('average energy (dB)'), ylabel('probability') ; 
title('GaussPDF for (en LF filt phAA)'); axis([-80 40 0 0.05]); grid on;
subplot(2,4,6);
plot(enLF_andHF_filt_phAA(:,2),p_HF_filt_phAA, '.');
xlabel('average energy (dB)'), ylabel('probability') ; 
title('GaussPDF for (en HF filt phAA)'); axis([-80 40 0 0.05]); grid on;
subplot(2,4,7);
plot(enLF_andHF_filt_phS(:,1),p_LF_filt_phS, '.');
xlabel('average energy (dB)'), ylabel('probability') ; 
title('GaussPDF for (en LF filt phS)'); axis([-80 40 0 0.05]); grid on;
subplot(2,4,8);
plot(enLF_andHF_filt_phS(:,2),p_HF_filt_phS, '.');
xlabel('average energy (dB)'), ylabel('probability') ; 
title('GaussPDF for (en HF filt phS)'); axis([-80 40 0 0.05]); grid on;


% function of Gaussian PDF, from tutorial
function [g_pdf]=gaussPDF(x,mu,var)

k= size(x);
g_pdf= ones(k) * (1/sqrt(2*pi*var));
Mu= ones(k) * mu;
Var= ones(k) * var;
g_pdf= g_pdf.* (exp(-((x-Mu).* (x-Mu))./(2*Var)));
end

