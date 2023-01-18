function [] = Main2()
%For running the second part of the exercise

%M-PAM

m_4=4;
m_8=8;

Tsymbol=4*10^(-6);
Tc=0.4*10^(-6);
Tsample=0.1*10^(-6);
fc=1/Tc;

bits = 6*10^5;

xsig = randsrc(1,bits,[0 1]);

[BER_4,SER_4] = MPAM(m_4,Tsymbol,bits,Tsample,fc,'',xsig);

[BER_8_,SER_8]=MPAM(m_8,Tsymbol,bits,Tsample,fc,'',xsig);


%{
m=8;

[BER_8,SER_8] = MPAM(m,Tsymbol,bits,Tsample,fc)

BER_8

SER_8

%}

SNR = 0:2:20;

%Plots of BER and SER for 4-PAM
figure(1)
semilogy(SNR,SER_4,'--',SNR,SER_8);
legend('SER 4-PAM','SER 8-PAM');
title(sprintf('Plot SER 4-PAM vs SER 8-PAM'));
xlabel('Signal Noise Ratio (SNR)');
ylabel('Symbol Error Rate (SER)');


%{
figure(2)
semilogy(0:2:20,SER_4,'-s');
title(sprintf('Plot for SER-SNR %d-PAM',m));
xlabel('Signal Noise Ratio (SNR)');
ylabel('Symbol Error Rate (SER)');
hold on;

%}

end

