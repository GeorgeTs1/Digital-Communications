function [xsig] = ADM_Decoder(ysig,d)

    sz = size(ysig,2);
        
    eq = zeros([1,sz]);
    
    xq = zeros([1,sz]);
    
    d(1) = eps;
    
    for i=1:sz
            
            if i==1

                eq(i) = ysig(i) * d(1);

                xq(i) = eq(i);

                xsig(i) = xq(i);

            else

           
                if ysig(i) == ysig(i-1)
                    
                    d(i) = d(i-1)*1.5;
                    
                else
                    
                    d(i) = d(i-1)/1.5;
                end
                
                
                eq(i) = ysig(i) * d(i);
                
                xq(i) = eq(i) + xq(i-1);
                
                xsig(i) = xq(i);   

            end

            
    end
    
        
end
    
    
    

  

