function [xi, Jmi] = single_segment_map(L, qi)

    if abs(qi) < 1e-6
      half_chord = L/2;
      d_half_chord = 0;  
    else
      half_chord = (L*sin(qi/2))/qi;
      d_half_chord = L*( qi*cos(qi/2) - 2*sin(qi/2) ) / (2*qi^2);
    end

    xi = [ qi/2;
           half_chord;
           half_chord;
           qi/2 ];

    Jmi = [ 1/2;
            d_half_chord;
            d_half_chord;
            1/2 ];
end
