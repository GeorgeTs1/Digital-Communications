function [BER,SER] = MPAM(m,Tsymbol,bits,Tsample,fc,encoding,xsig)


symbols = log2(m); % Symbols of M-PAM

%xsig = randsrc(1,bits,[0 1]) % Creation of equiprobable 0 and 1

x_size = size(xsig,2); %Size of signal (columns)

%Each row has a log2(m) binary digits will then be input to mapper
for i=1:symbols
    
    source(:,i) = xsig(i:symbols:x_size);

end

cnt = 1; % For BER and SER purposes

t = 0:Tsample:Tsymbol; % Time-Period of carrier signal

Es = 1; %Energy per Symbol

gt = sqrt(2/Tsymbol); % Rectangualar Pulse


decade_source = bi2de(source,'left-msb'); %Conversion to decade source w.r.t left-msb

if strcmp(encoding,'gray') == 1

    decade_source = bin2gray(decade_source,'pam',m);


end


signal = {}; % Cell array for storing the sm waveforms

r = {}; % Cell array for storing the output of AWGN channel waveforms

est_symbol_snr = {}; % Cell array for estimated symbol per SNR

est_binary_snr = {}; % Cell array for estimated binary sequence per SNR


%For each SNR in range of [0,20]
for SQNR = 0:2:20

    input = decade_source'; %For convenient purposes

    %Modulation
    for i=1:length(input)
                
        amplitude = (2*(input(i)+1)-(m+1)); %Amplitude of each symbol
            
       %Signal Function storing each waveform of each symbol sm
      
       signal{i} = amplitude*gt*cos(2*pi*fc*t);


       %Noise Calculation based on the formula SNR = 10*log10(Pm/Pn)
       power=sum(signal{i}.^2)/(length(signal{i}));

       noise_var = 10*power/(10^(SQNR/10));
       
       noise= sqrt(noise_var)*randn(1,length(signal{i}));

       % Output of AWGN channel for each symbol waveform sm
       r{i} = noise + signal{i};


    end 


     %Demodulation
    
     %Detector
        for k=1:length(input)
      
            f{k}  = r{k}.*gt.*cos(2*pi*fc*t)*Tsample;

          
            fsum(k) = sum(f{k});
        
        end
     
      %Demapper
    
        for l=1:length(input)

           minimum = inf;

            for j=1:m

                A = 2*j-(m+1);

                dist = abs(A-fsum(l));
               
                if (dist<minimum) 

                    minimum = dist;
                    right = j; 

                end
                

            end

             est_symbol(l) = right-1;

        end

        est_symbol_snr{cnt} = est_symbol;

        if strcmp(encoding,'gray')

               est_symbol = gray2bin(est_symbol,'pam',m);

        end

        est_binary =  de2bi(est_symbol,'left-msb')';

  
        est_binary_snr{cnt} = est_binary(:);



        cnt = cnt+1;

        
end


%BER and SER calculations

%BER
for i=1:length(est_binary_snr)

   [~,BER(i)]=biterr(est_binary_snr{i},xsig');
 
end

%SER

for j=1:length(est_symbol_snr)

   [~,SER(j)] = symerr(est_symbol_snr{j},decade_source');
   

end


end
    