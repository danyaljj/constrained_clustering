function [v,vi,nvi,greedy,greedy_11,hc,hk,h_ck,h_kc] =  calculate_measures (matrix)

%input: matrix - a co-occurrence matrix between induced and gold clusters
%lines are for induced clusters, columns are for gold clusters
  
% output: v,vi,nvi,greedy many-to-one (in the greedy variable) and greedy
% one-to-one (in the greedy_11 variable).
% hc and hk - unconditional entropies of gold (hc) and induced (hk) clusterings
% h_ck and h_kc: conditional entropies: h_ck - gold given induced. h_kc:
% induced given gold

  

         
  
   [v,hc,hk,h_ck,h_kc] = calculate_v_measure (matrix'); 
   
   vi = calculate_vi_measure (matrix');
  
   nvi = calculate_nvi_measure (matrix');
    
   greedy = calculate_greedy (matrix');

   greedy_11  = calculate_hungarian_munkers (matrix);
 
