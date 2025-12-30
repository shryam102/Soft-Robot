function [xi, Jmi] = single_map(qi)
phi   = qi(1);
theta = qi(2);

[li, dl_theta, dl_S] = Link_func(qi);

xi = [   phi;
      theta/2;
          li;
      theta/2;
        -phi];

Jmi = [1      0       0;
       0     1/2      0;
       0   dl_theta  dl_S;
       0     1/2      0;
      -1      0       0];

end