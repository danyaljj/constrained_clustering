function [v,hc,hk,h_ck,h_kc] = calculate_v_measure (matrix)
  
  %I assume that rows represent gold clusters and columns represent
  %induced clusters
  
  
 rows_number    = size (matrix,1);
 columns_number = size (matrix,2);


 rows_sum    = zeros (rows_number,1);
 columns_sum = zeros (columns_number,1);
 N           = sum(sum(matrix));
 
 
 for i=1:rows_number
   rows_sum(i) = sum(matrix(i,:));
 end
   
   
 for i=1:columns_number
   columns_sum(i) =  sum(matrix(:,i));
 end  
   
 
  h_ck = 0; 
  h_kc = 0;
  hk   = 0;
  hc   = 0;
   
   for i=1:rows_number
     for j=1:columns_number
       
       if (matrix (i,j) == 0)
         continue;
         
       end
       
       h_ck = h_ck + (matrix (i,j)/N) * log (matrix(i,j) / columns_sum(j));
       h_kc = h_kc + (matrix (i,j)/N) * log (matrix(i,j) / rows_sum(i));
      
     end
   end
   
   h_ck = - h_ck;
   h_kc = - h_kc;
   
   
   
      
   
   for j = 1:rows_number
     if (rows_sum(j) > 0)
         hc = hc + (rows_sum(j)/N) * log (rows_sum(j) / N);
     end
   end
   
   
   
   for j = 1:columns_number
      if (columns_sum (j) > 0) 
        hk = hk + (columns_sum(j)/N) * log (columns_sum (j)/N);
      end
   end
     
     
     hc = -hc;
     hk = -hk;


   h = 0;
   c = 0;
   
   hc;
   hk;
   h_ck;
   h_kc;
   
   if (hc == 0)
     h = 1;
   else
     h = 1 - (h_ck/hc);
   end
   
     
     if (hk == 0)
       c = 1;
     else
       c = 1 - (h_kc/hk);
     end 
     
     h;
     c;
     
   if ( (h == 0) | (c == 0) )
        v = 0;
   else    
     v = 2*h*c/(h+c);
     
   end
   
   