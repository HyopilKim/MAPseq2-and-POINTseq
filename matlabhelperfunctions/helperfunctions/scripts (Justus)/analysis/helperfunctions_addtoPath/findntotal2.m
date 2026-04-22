%% function to derive N_total,as described in the MAPseq methods section
function [N_t,matrix_binary]=findntotal2(matrix,cutoff)

%binarize input
matrix_binary=matrix>cutoff;
matrix_binary=matrix_binary(sum(matrix_binary,2)~=0,:);

%find counts
N_a=sum(matrix_binary,1);

[N_obs, k] = size(matrix_binary);
   N_a = sum(matrix_binary, 1);  % count 1s in each dimension
   
   % Initialize polynomial coefficients array
   coeffs = zeros(1, k+1);
   
   % First coefficient (highest degree)
   coeffs(1) = N_obs - sum(N_a);
   
   % Calculate other coefficients using nchoosek
   for i = 2:k
       % (-1)^(i-1) gives alternating signs
       coeffs(i) = (-1)^(i-1) * sum(prod(nchoosek(N_a,i),2));
   end
   
   % Last coefficient (constant term)
   coeffs(end) = (-1)^k * prod(N_a);
   
   % Find roots
   r = roots(coeffs);
   
   % Select real and positive roots
   real_positive_roots = r(imag(r)==0 & real(r)>0);
   
   % Take the smallest positive real root as N_total
   N_total = min(real_positive_roots);