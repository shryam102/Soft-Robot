function [xi, Jm] = build_xi_Jm(L, q)
  n   = numel(q);
  xi  = zeros(4*n,1);
  Jm  = zeros(4*n,n);
  for i = 1:n
    idx = (i-1)*4 + (1:4);  
    [xi_block, Jmi] = single_segment_map(L(i), q(i));
    xi(idx)        = xi_block;     
    Jm(idx, i)     = Jmi;          
  end
end
