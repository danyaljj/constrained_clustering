function score = calculate_hungarian_munkers(M)


    size(M);
  
    k = max(max(M));
    
    k = k + 100;
    
    W = k - M;
    
    
  [A,cost] = assignmentoptimal(W);

  A
  
  size(A);
  
  l = size(M,1);
  
  sum1 = 0;
  
  for i =1:1:l
    
    sum1 = sum1 + M(i,A(i));
    
  end
    
  sum1
  
  score = sum1/sum(sum(M))
