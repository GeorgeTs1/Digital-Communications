function [] = Main()


b=1;

a1 = 0.9; % coefficient for  AR1(1)

a2 = 0.95; % coefficient for  AR2(1)




a1_=[1,-a1]; %AR1(1)

a2_=[1 -a2]; %AR2(1)

L = 10000; %Size of white noise

x = randn(1,L); % White noise transpose

y1 = filter(b,a1_,x); % AR1(1)

y2 = filter(b,a2_,x); % AR2(1)

k=1; % For Illustration purposes

%Lloyd Max for AR1(1) AR2(1)
for i=1:3

 N = 2^i;   

figure(k);
[xq1,centers1,D1] = LloydMax(y1,N,min(y1),max(y1),1);
figure(k+1);
[xq2,centers2,D2] = LloydMax(y2,N,min(y2),max(y2),2);

%Entropy Calculation Process
p1 = hist(xq1,N);
p2 = hist(xq2,N);

p1(p1==0) = []; % Remove zero entries
p1 = p1 ./ numel(xq1);

H1(i) = -sum(p1.*log(p1));

p2(p2==0) = []; %Remove Zero Entries
p2 = p2 ./ numel(xq2);

H2(i) = -sum(p2.*log(p2));

k = k+2;

end


H1

H2

%Image Processing
img_src=load('cameraman.mat'); 
 
img = img_src.i; %Loading the matrix of pixels
 
figure(k),imshow(uint8(img)); %Image Projection
 
x = img(:); %One column vector of values
 
x = (x-128)/128; %Scaling from range [0 255] to [-1 1]
 
k = k+1;
 
y = cell(1,2); %For storing the images

%Lloyd-Max for Image 256*256
for i=1:2
    
    N = 2^i; 
    
    figure(k);
    
    [xq,centers,D] = LloydMax(x',N,min(x'),max(x'),3);
    
    %Entropy Calculation Process
    p = hist(xq,N);
    
    p(p==0) = []; % Remove zero entries
    p = p ./ numel(xq);

    H_img(i) = -sum(p.*log(p)); %Entropy of Image

    xq = xq*128 + 128; %Scaling up to [0 255] from [-1 1]
    
    y{1,i} = reshape(xq',256,256); %Reshape vector to matrix
    
    k = k+1;
   
    
end


H_img

figure(k),

subplot(1,2,1),imshow(uint8(y{1,1}));
title('PCM 2 bits');
subplot(1,2,2),imshow(uint8(y{1,2}));
title('PCM 4 bits');

k=1;

%ADM 
for i=1:3
    
M = 2^i;

%Oversampling by a factor of M
y1_sampled = interp(y1,M);
y2_sampled = interp(y2,M);

[y1_,d1] = ADM_coder(y1_sampled); %ADM_coder for AR1(1)
x1  = ADM_Decoder(y1_,d1); %ADM_decoder for AR1(1)


[y2_,d2] = ADM_coder(y2_sampled); %ADM_coder for AR2(1)
x2 = ADM_Decoder(y2_,d2); %ADM_decoder for AR2(1)

% SNR1 for AR1(1) and SNR2 for AR2(1) 
SNR1(k) = 10*log10(mean(y1_sampled.^2)/mean((y1_sampled-x1).^2)); %SNR for source AR1(1)
SNR2(k) =  10*log10(mean(y2_sampled.^2)/mean((y2_sampled-x2).^2)); % SNR for source AR2(1)

k = k +1;

end

SNR1

SNR2





end


