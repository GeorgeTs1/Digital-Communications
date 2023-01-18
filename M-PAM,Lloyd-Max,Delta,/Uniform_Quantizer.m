function [centers] = Uniform_Quantizer(x,N,min_val,max_val)

L = 2^N; %Quantization Levels

delta = (max_val-min_val)/L; %Step size 

centers = zeros([1,L+2]);

centers(1) = min_val;

centers(2) = min_val+delta/2;

%Calculating the centers of Uniform Quantizer
for i=3:L+2
    centers(i) = centers(i-1) + delta; 
end

    centers(i) = max_val;

end