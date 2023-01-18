function [] = Main()
%Main for M-PSK and Ploting purposes


bits = 6*10^5;

src = randsrc(1,bits,[0 1]);

[BER_4,SER_4,BER_h1_4,SER_h1_4,BER_h2_4,SER_h2_4]=MPSK(4,src);
[BER_8,SER_8,BER_h1_8,SER_h1_8,BER_h2_8,SER_h2_8]=MPSK(8,src);

snr = 0:2:30;

%Plots of BER and SER
figure(1)
semilogy(snr,BER_4,snr,BER_h1_4,snr,BER_h2_4,snr,BER_8,snr,BER_h1_8,snr,BER_h2_8);

title('Plot for BER-SNR PSK');
xlabel('Signal Noise Ratio (SNR)');
ylabel('Bit Error Rate (BER)');
legend('ideal(4-PSK)','h1(4-PSK)','h2(4-PSK)','ideal(8-PSK)','h1(8-PSK)','h2(8-PSK)')

figure(2)

semilogy(snr,SER_4,snr,SER_h1_4,snr,SER_h2_4,snr,SER_8,snr,SER_h1_8,snr,SER_h2_8);

title('Plot for SER-SNR PSK');
xlabel('Signal Noise Ratio (SNR)');
ylabel('Symbol Error Rate (SER)');
legend('ideal(4-PSK)','h1(4-PSK)','h2(4-PSK)','ideal(8-PSK)','h1(8-PSK)','h2(8-PSK)')


%Constellation Plots

ref_4 = [0:3];
ref_8 = [0:7];

ref_gray_4 = bin2gray(ref_4,'psk',4);
ref_gray_8 = bin2gray(ref_8,'psk',8);
   
ipBin_4 = dec2bin(ref_gray_4.');
ipBin_8 = dec2bin(ref_gray_8.');

ipPhase_4 = ref_4*2*pi/4;
ipPhase_8 = ref_8*2*pi/8;

mod_4 = exp(1j*ipPhase_4);
mod_8 = exp(1j*ipPhase_8);

figure(3)

scatter(real(mod_4),imag(mod_4));

text(real(mod_4)-0.1,imag(mod_4)-0.1,ipBin_4);

axis([-1.25 1.25 -1.25 1.25]);

grid on

figure(4)

scatter(real(mod_8),imag(mod_8));

text(real(mod_8)-0.1,imag(mod_8)-0.1,ipBin_8);

axis([-1.25 1.25 -1.25 1.25]);

grid on

end