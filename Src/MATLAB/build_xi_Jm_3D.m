function [xi, Jm] = build_xi_Jm_3D(L, q)
  n   = numel(q)/2;
  xi  = zeros(10*n,1);
  Jm  = zeros(10*n,2*n);
  for i = 1:n
    idx = (i-1)*10 + (1:10);  
    [xi_block, Jmi] = single_segment_map_3D(L(i), q(2*i-1:2*i));
    xi(idx)        = xi_block;     
    Jm(idx, 2*i-1:2*i)     = Jmi;          
  end
end
