function [ysig,d] = ADM_coder(xsig)

  

if size(xsig,1) ~= 1
      error('Please enter x with dimensions [1,length(xsig)]');
end



sz = size(xsig,2); % size of signal

bn = zeros([1,sz]); % Initialize of b

eq = zeros([1,sz]);

xq = zeros([1,sz]);

e = zeros([1,sz]);

ysig = zeros([1,sz]); %Initialize output

d = zeros([1,sz]); % Initialize of d

d(1) = eps;        %Smallest value possible

i=1;


while i<=sz
    
    if i==1
       
        e(1) = xsig(1);
        
    else
       
        e(i) = xsig(i) - xq(i-1);
        
        
    end
    
     if(e(i)<0)

            b(i) = -1;
           
            if i==1 % Only for the first time
               
               eq(i) = b(i) * d(i);
               
               xq(i) = eq(i);

            else
                
              if b(i)==b(i-1) % Update the step 
                
                   d(i) = d(i-1) * 1.5;
                   
                   eq(i) = b(i) * d(i);

              else
                    
                    d(i) = d(i-1)/1.5;
                    
                    eq(i) = b(i) * d(i); 

              end

              
               xq(i) = eq(i) + xq(i-1);

            end

     else

         b(i) = 1;
            
         if i==1
            
               eq(i) = b(i) * d(i);
               
               xq(i) = eq(i);

         else
                
              if b(i)==b(i-1)
                
                   d(i) = d(i-1) * 1.5;
                  
                   eq(i) = b(i) * d(i);

              else
                    d(i) = d(i-1)/1.5; 
                  
                    eq(i) = b(i) * d(i); 

              end

              
               xq(i) = eq(i) + xq(i-1);

            end

     end

     ysig(i) = b(i);

     i = i +1;

end



end