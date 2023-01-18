function [xq,centers,D] = LloydMax(x,N,min_val,max_val,sel)

%If sel==1 we select the AR1(1) source else AR2(1)

if sel==1
    
    source = 'AR1(1)';
    
elseif sel==2
    
    source = 'AR2(1)';
    
else
    
    source = 'Image';
    
end



if size(x,1) ~= 1
    error('Please enter x with dimensions [1,length(x)]');
end


L = 2^N; % Quantization Levels

centers = Uniform_Quantizer(x,N,min_val,max_val); %Centers of uniform quantizer

iters = 1;

D = []; %Distortion Vector Initialization

SQNR = []; %SQNR Vector for each iteration of the Algorithm

%Main Loop of the algorithm
while true
    
    total = 0;
    
    T = []; % Quantization Zones Vector
    
    cond_mean = zeros([1,length(centers)]); % Conditional Mean Vector
    
    counter = zeros([1,length(centers)]); % Counter of each x(i) in range Tk < x(i) < Tk+1
    
    T(1) = min_val; %First value is always the minimum of signal
    %Calculate the Quantization Zones
    for i=2:(length(centers)-2)
        
        T(i) = (centers(i) + centers(i+1))/2;
    end
    
    T(i+1) = max_val; %Last value is always the maximum of signal
    
    
    %Quantization Process
    for i=1:length(x)
        
        for j=1:(length(T)-1)
            
            if x(i)>T(j) && x(i)<=T(j+1)
                
                xq(i) =  centers(j+1);
                
                %Mean distance from the center
                total = total + abs(centers(j+1) - x(i)).^2;
                
                %Conditional Mean of zones
                cond_mean(j) = cond_mean(j) + x(i);
                
                counter(j) = counter(j) + 1;
                
            end
            
        end
        
        if x(i) == T(1)
            
            xq(i) = centers(2);
            
            total = total + abs(centers(2) - x(i)).^2;
            
            cond_mean(1) = cond_mean(1) + x(i);
            
            counter(1) = counter(1) + 1;
            
        end
    end
    
    
    avg_distortion = total/length(x);
    
    D(iters) = avg_distortion;
    
    
    SQNR(iters) = 10*log10(mean(x.^2)/D(iters)); %SQNR Each Iteration
    
    
    if SQNR(iters) == inf % Meaning no distortion
        
        SQNR(iters) = realmax;
        
    end
    
    
    
    %Update the centroids
    for j=2:(length(centers)-1)
        
        if counter(j-1) ~= 0
            centers(j) = cond_mean(j-1) / counter(j-1);
        end
    end
    
    % Termination criterion of loop
    if length(D) ~= 1
        if abs(D(iters) - D(iters-1)) < eps
            
            break;
            
        end
        
    end
    
    iters = iters + 1;
    
    
    
    
end





centers(1)= [];
centers(end) = [];

iters = [1:1:iters];



%SQNR Plot for each iteration
subplot(2,2,[1,2])
plot(iters,SQNR);
title(sprintf('SQNR-Iterations plot for %d bits quantizer %s',N,source));
xlabel('Iterations');
ylabel('SQNR(DB)');


%Original signal plot
subplot(2,2,3)
plot(x)
title(sprintf('Original waveform %d bits quantizer %s',N,source));


%Output of quantized signal
subplot(2,2,4)
plot(xq);
title(sprintf('Output waveform %d bits quantizer %s',N,source));




end


