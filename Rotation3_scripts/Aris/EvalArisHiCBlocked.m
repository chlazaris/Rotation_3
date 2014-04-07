function [E] = EvalArisHiCBlocked(f1,f2,fchr)

  % read input matrices
  C1 = sparse(load(f1)); 
  C2 = sparse(load(f2)); 
  n = 200;
  C1 = C1(1:n, 1:n);
  C2 = C2(1:n, 1:n);

  % Set constants
  gamma = 0.1;
  epsilon = 1e-6;
  maxit = 100;
  theta = 0.1;
  cp_map = [1 n+1];
  normal = 0;
  len = 0;
  gc = 0;
  map = 0;



  % Calculate Spearman correlation between P1 and P2
  
  L = 0:0.1:10;
  T = 0:0.1:10;

  E = zeros(length(L),length(T));

  for k_map=10;
    for i_lambda = 1:length(L);
      lambda = L(i_lambda);
      for i_tau = 1:length(T);
        tau = T(i_tau);
        fprintf('%d,%d\n',lambda,tau);

        P1 = ArisHiCBlocked(C1, tau,  gamma,  lambda,  epsilon,  maxit,  theta,  cp_map,  k_map,  normal,  len,  gc,  map);
        P2 = ArisHiCBlocked(C2, tau,  gamma,  lambda,  epsilon,  maxit,  theta,  cp_map,  k_map,  normal,  len,  gc,  map);
        
        % Turn the matrices to vectors
        E(i_lambda,i_tau) = corr(P1(:), P2(:), 'type', 'Pearson');

       end
     end
   end

