%
%function [P]= ArisHiCBlocked(C,...
%                         tau,...
%                         gamma,...
%                         epsilon,...
%                         maxit,...
%                         theta,...
%                         cp_map,...
%                         k_map,...
%                         normal,...
%                         len,...
%                         gc,...
%                         map)
% HiC Solver using Incremental Proximal Descent
%
% Parameters:
% C --> observed contact matrix
% tau --> penalty for nuclear norm
% gamma --> penalty for regression coefficients
% lambda --> weight for L1(P) in objective function
% epsilon --> tolerance for convergence
% maxit --> maximum number of iterations 
% theta --> parameter for gradient descent direction
% cp_map --> cumulative partition map that determines the diagonal blocks
% k_map --> map for the rank-k approximation for each partition.
% normal --> Should we normalize the diagonal blocks or not
% len -> Vector containing the length of each bin (nx1)
% gc --> GC content of each bin (nx1)
% map -> map content of each bin (nx1)
%
% Return:
% P (sparse) estimated contact matrix
%

function [P]= ArisHiCBlocked(C,...
                         tau,...
                         gamma,...
                         lambda,...
                         epsilon,...
                         maxit,...
                         theta,...
                         cp_map,...
                         k_map,...
                         normal,...
                         len,...
                         gc,...
                         map)

% STDOUT is often not defined in MATLAB
if (0 == exist('stdout','var')) stdout=1; end

% Make sure we are dealing with column vectors
if (1 == size(cp_map,1)) cp_map = cp_map'; end
if (1 == size(k_map,1)) k_map = k_map'; end

% Get the sizes
n = size(C,1);
nparts = (size(cp_map,1)-1);
 
% Error 1 --- The partition and rank-k maps should have the same size
if (size(cp_map,1) ~= (1+size(k_map,1))) 
  error('Partition and rank-k map size mismatch (%d,%d)\n', ...
         nparts, size(k_map,1));
end

% Error 2 --- The rank-k approximation requested should be < partition size
for i=1:nparts
  partition_size = (cp_map(i+1)-cp_map(i)+1);
  if (k_map(i) >= partition_size) 
    error('Invalid low-rank approximation requested (%d,%d)\n',...
          partition_size, k_map(i));
  end
end

% Error 3 --- The sum of all the partitions should add up to n
if ((n+1) ~= cp_map(nparts+1))
  error('Paritition map does not add up to size of the matrix\n');
end

% Print out all the function parameters
fprintf (stdout, '-------------------------------\n');
fprintf (stdout, 'HicBlocked function parameters\n');
fprintf (stdout, 'tau = %f\n', tau);
fprintf (stdout, 'gamma = %f\n', gamma);
fprintf (stdout, 'lambda = %f\n', lambda);
fprintf (stdout, 'epsilon = %f\n', epsilon);
fprintf (stdout, 'theta = %f\n', theta);
fprintf (stdout, 'maxit = %d\n', maxit);
fprintf (stdout, 'normalize C = %d\n', normal);
fprintf (stdout, 'cp_map %f\n', cp_map');
fprintf (stdout, 'k_map %f\n', k_map');
fprintf (stdout, '-------------------------------\n');

% First normalize if requested 
if (normal) C = HiCHarvard (C, cp_map, len, gc, map); end

% Decompose C into diagonal and off-diagonal blocks. 
[C_diag,C_off_diag] = SplitSymMatrix(C,cp_map);

% Our initial guess P=C
P_diag = C_diag;
P_off_diag = C_off_diag;

% Incremental Proximal Descent
for ite=1:maxit 

  l1_of_S = 0;
  l1_of_P = 0;
  fro_of_P_minus_C = 0;

  for i=1:nparts

    % Gradient Descent --- P=P-(2*theta*(P-C)); 
    P_diag{i} = P_diag{i} - (2*theta*(P_diag{i}-C_diag{i}));
    P_off_diag{i} = P_off_diag{i} - (2*theta*(P_off_diag{i}-C_off_diag{i}));

    % Low rank approximate the diagonal part
    if tau>0
        [U,S,V] = irlba(P_diag{i}, struct('K',k_map(i)));
        S = (S-tau) .* ((S-tau)>0);
        UTU = ((U*S*U')-gamma);
        P_diag{i} = sparse(UTU.*(UTU>0));
    else
        P_diag{i} = (P_diag{i}-gamma);
        P_diag{i} = sparse(P_diag{i}.*(P_diag{i}>0));
    end

    % Accumulate terms for the objective function from diagonal part
    if tau>0 l1_of_S = l1_of_S + sum(sum(S)); end
    l1_of_P = l1_of_P + sum(sum(P_diag{i}));
    fro_of_P_minus_C = fro_of_P_minus_C + ...
                       sum(sum((P_diag{i}-C_diag{i}).^2));

    % Threshold the off-diagonal part
    P_off_diag{i} = P_off_diag{i}-gamma;
    P_off_diag{i} = P_off_diag{i}.*(P_off_diag{i}>0);

    % Accumulate terms for the objective function from diagonal part
    l1_of_P = l1_of_P + 2*sum(sum(P_off_diag{i}));
    fro_of_P_minus_C = fro_of_P_minus_C + ...
                       2*sum(sum((P_off_diag{i}-C_off_diag{i}).^2));

  end

  % Compute the objective function
  L = fro_of_P_minus_C + lambda*l1_of_P + tau*l1_of_S;

  fprintf (stdout, 'Iteration = %d; fro(P-C)^2 = %f; L1(P) = %f; objective = %f\n', ite, fro_of_P_minus_C, l1_of_P, L);

  % Check for stopping condition
  if ((ite > 1) && ((abs(L-prevL)/abs(prevL)) < epsilon)) break; end

  % Save this iteration's objective function
  prevL=L;
  
end

P = JoinMatrix(P_diag, P_off_diag, cp_map, n);
