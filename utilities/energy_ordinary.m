function loop_energy = energy_ordinary()

% ENERGY_ORDINARY  Thin-plate energy for an ordinary Loop surface patch.

loop_energy = (1/486) * ...
      [   5    3   -1    2   -9  -15    1   -6    3    9    1    2 ;
          0   21    3   -6  -33  -33   -6    3    6    3    9    9 ;
          0    0    5    1  -15   -9    2    9    3   -6    2    1 ;
          0    0    0    5   -9    3    2    3  -15    9   -1    1 ;
          0    0    0    0   57   -6    3  -33   -6    6  -15    3 ;
          0    0    0    0    0   57   -9    6   -6  -33    3  -15 ;
          0    0    0    0    0    0    5    9  -15    3    1   -1 ;
          0    0    0    0    0    0    0   21  -33    3    3   -6 ;
          0    0    0    0    0    0    0    0   57  -33   -9   -9 ;
          0    0    0    0    0    0    0    0    0   21   -6    3 ;
          0    0    0    0    0    0    0    0    0    0    5    2 ;
          0    0    0    0    0    0    0    0    0    0    0    5 ];

loop_energy = loop_energy + loop_energy';

end

