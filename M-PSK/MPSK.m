function [BER,SER,BER_h1,SER_h1,BER_h2,SER_h2] = MPSK(M,src)
%M-PSK Modulation and Demodulation 


warning('off');

symbols = log2(M); % Symbols 

roll_off = 0.3; % For square root raised cosine
sps = 4; %Samples per symbol
span = 6; %Truncation to 6 symbols

Es = 1; %Energy per Symbol

Tsymbol = 20;

Tsample = 4;

Tc = 24;

fc = 1/Tc;

t = 0:Tsample:Tsymbol-Tsample;

h_f = rcosdesign(roll_off,span,sps,'sqrt'); %Square root raised cosine
 
N = length(src)/symbols; % Number of symbols of  binary sequence

%Channels (non ideal)
h1 = [0.04 -0.05 0.07 -0.21	-0.5 0.72 0.36 0 0.21 0.03 0.07];

h2 = [0.227 0.460 0.688 0.460 0.227];

h1 = upsample(h1,4); %Upsample the impulse response of Channel 1(non ideal) 

h2 = upsample(h2,4); %Upsample the impulse response of Channel 2(non ideal)


%Grouping the symbols in log2(M) teams 

   x_mod = reshape(src,symbols,N)';

   x_mod_dec = bi2de(x_mod,'left-msb');  %Binary to decimal conversion
    
   %Gray Encoding
   x_gray = bin2gray(x_mod_dec,'psk',M); %Gray encoding
    
       
cnt = 1; % For BER and SER indexing

%Mod Symbols Waveforms 
x_mod = cell(1,length(src)/symbols);
x_mod_h1 = cell(1,length(src)/symbols);
x_mod_h2 = cell(1,length(src)/symbols);

% Component Vectors for Receiver
r1 = cell(1,length(src)/symbols);
r2 = cell(1,length(src)/symbols);

r1_h1 = cell(1,length(src)/symbols);
r2_h1 = cell(1,length(src)/symbols);

r1_h2 = cell(1,length(src)/symbols);
r2_h2 = cell(1,length(src)/symbols);

%Estimate Symbols Vectors
est_symbol = zeros(1,length(src)/symbols);
est_symbol_h1 = zeros(1,length(src)/symbols);
est_symbol_h2 = zeros(1,length(src)/symbols);

%BER SER vectors initialization
BER = zeros(1,16); %Number of 0:2:30
BER_h1 = zeros(1,16); %Number of 0:2:30
BER_h2 = zeros(1,16); %Number of 0:2:30

SER = zeros(1,16); %Number of 0:2:30
SER_h1 = zeros(1,16); %Number of 0:2:30
SER_h2 = zeros(1,16); %Number of 0:2:30

for SNR = 0:2:30
       
    for k = 1:length(x_gray)
       
        phi = 2*pi*x_gray(k)/M;

        Amc = cos(phi);

        Ams = sin(phi);

        sm = sqrt(2/Tsymbol)*Amc*cos(2*pi*fc*t) - sqrt(2/Tsymbol)*Ams*sin(2*pi*fc*t);
 
        sm_up = upsample(sm,4);

        sm_filt = conv(sm_up,h_f,'same');

       %Noise calculation

       power=sum(sm_filt.^2)/(length(sm_filt));
            
       s=10*power/(10^(SNR/10));
       
       noise = sqrt(s)*randn(1,length(sm_filt));
       
       x_mod{k} =  sm_filt + noise;

       x_mod_h1{k} = conv(sm_filt,h1,'same') + noise;

       x_mod_h2{k} = conv(sm_filt,h2,'same') + noise;

    end

         %Detector
    for l=1:length(x_gray)

           %Filtering and downsampling
            tmp_x = conv(x_mod{l},h_f,'same');
            tmp_xx = downsample(tmp_x,4);

            tmp_x_h1 = conv(x_mod_h1{l},h_f,'same');
            tmp_xx_h1 = downsample(tmp_x_h1,4);

            tmp_x_h2 = conv(x_mod_h2{l},h_f,'same');
            tmp_xx_h2 = downsample(tmp_x_h2,4);

           %Mapping to the orthonormal basis functions of M-PSK 
            I = sqrt(2/Tsymbol) .* tmp_xx .* cos(2*pi*fc*t)*Tsample;
            Q = -sqrt(2/Tsymbol) .* tmp_xx .* sin(2*pi*fc*t)*Tsample;

            I_h1 = sqrt(2/Tsymbol) .* tmp_xx_h1 .* cos(2*pi*fc*t)*Tsample;
            Q_h1 = -sqrt(2/Tsymbol) .* tmp_xx_h1 .* sin(2*pi*fc*t)*Tsample;

            I_h2 = sqrt(2/Tsymbol) .* tmp_xx_h2 .* cos(2*pi*fc*t)*Tsample;
            Q_h2 = -sqrt(2/Tsymbol) .* tmp_xx_h2 .* sin(2*pi*fc*t)*Tsample;

            %Integration and production of components
            r1{l} = trapz(t,I);
            r2{l} = trapz(t,Q);

            r1_h1{l} = trapz(t,I_h1);
            r2_h1{l} = trapz(t,Q_h1);

            r1_h2{l} = trapz(t,I_h2);
            r2_h2{l} = trapz(t,Q_h2);
        
    end

   
        %Receiver 
    
    for m=1:length(x_gray)
      
       min_ind = inf;
       min_ind_h1 = inf;
       min_ind_h2 = inf;
   
       for j=1:M

               phase = 2*pi*(j-1)/M;

               dist = (r1{m}-cos(phase)).^2 + (r2{m}-sin(phase)).^2;
               dist_h1 = (r1_h1{m}-cos(phase)).^2 + (r2_h1{m}-sin(phase)).^2;
               dist_h2 = (r1_h2{m}-cos(phase)).^2 + (r2_h2{m}-sin(phase)).^2;

               if dist < min_ind

                    min_ind = dist;

                    right = j-1;

               end


               if dist_h1 < min_ind_h1

                    min_ind_h1 = dist_h1;

                    right_h1 = j-1;

               end

                if dist_h2 < min_ind_h2

                    min_ind_h2 = dist_h2;

                    right_h2 = j-1;

               end


       end
                 
            est_symbol(m) = right;
            est_symbol_h1(m) = right_h1;
            est_symbol_h2(m) = right_h2;

   end



    est_symbols = gray2bin(est_symbol,'psk',M);
    est_symbols_h1 = gray2bin(est_symbol_h1,'psk',M);
    est_symbols_h2 = gray2bin(est_symbol_h2,'psk',M);

   

    est_bin = de2bi(est_symbols,symbols,'left-msb')';
    est_bin_h1 = de2bi(est_symbols_h1,symbols,'left-msb')';
    est_bin_h2 = de2bi(est_symbols_h2,symbols,'left-msb')';

    est_bin = est_bin(:);
    est_bin_h1 = est_bin_h1(:);
    est_bin_h2 = est_bin_h2(:);


    %BER SER

    [~,BER(cnt)] = biterr(src,est_bin');
    [~,BER_h1(cnt)] = biterr(src,est_bin_h1');
    [~,BER_h2(cnt)] = biterr(src,est_bin_h2');

    [~,SER(cnt)] = symerr(x_mod_dec',est_symbols);
    [~,SER_h1(cnt)] = symerr(x_mod_dec',est_symbols_h1);
    [~,SER_h2(cnt)] = symerr(x_mod_dec',est_symbols_h2);


    cnt = cnt+1;
    
end 

end
