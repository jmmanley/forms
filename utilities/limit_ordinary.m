function [ limitpoint deriv secderiv ] = limit_ordinary( u, v )

% LIMIT_ORDINARY  Find limit point on a subdivision surface
%   where the required point lies in an ordinary surface patch (i.e. next
%   to three regular vertices)
%
%   See also LIMIT_EXTRAORDINARY

w = 1 - (u + v);

    % This matrix is from
    %   Evaluation of Piecewise Smooth Subdivision Surfaces
    %   Denis Zorin, Daniel Kristjansson
    %   The Visual Computer, 2002
    %
    % They say they got it from
    %   Analysis and Application of Subdivision Surfaces
    %   J. E. Schweitzer
    %   PhD thesis, University of Washington, Seattle, 1996
    tobezier = [ 2  2  0  2 12 2  0  2 2  0  0  0 ;
                 0  1  0  1 12 3  0  3 4  0  0  0 ;
                 1  3  0  0 12 4  0  1 3  0  0  0 ;
                 0  0  0  0 8  4  0  4 8  0  0  0 ;
                 0  1  0  0 10 6  0  1 6  0  0  0 ;
                 0  4  0  0 8  8  0  0 4  0  0  0 ;
                 0  0  0  0 4  3  0  3 12 1  1  0 ;
                 0  0  0  0 6  6  0  1 10 1  0  0 ;
                 0  1  0  0 6  10 0  0 6  1  0  0 ;
                 0  3  1  0 4  12 0  0 3  1  0  0 ;
                 0  0  0  0 2  2  0  2 12 2  2  2 ;
                 0  0  0  0 3  4  0  1 12 3  0  1 ;
                 0  0  0  0 4  8  0  0 8  4  0  0 ;
                 0  1  0  0 3  12 1  0 4  3  0  0 ;
                 0  2  2  0 2  12 2  0 2  2  0  0 ] ./ 24;

    rearrange = [ 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ;
                  0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 ;
                  0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 ;
                  0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 ;
                  0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 ;
                  0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 ;
                  0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 ;
                  0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 ;
                  0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 ;
                  0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 ;
                  0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 ;
                  0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 ;
                  0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 ;
                  0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 ;
                  0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 ];

    tobezier = rearrange * tobezier;

    decast  = [ v u 0 0 0 w 0 0 0 0 0 0 0 0 0 ;
                0 v u 0 0 0 w 0 0 0 0 0 0 0 0 ;
                0 0 v u 0 0 0 w 0 0 0 0 0 0 0 ;
                0 0 0 v u 0 0 0 w 0 0 0 0 0 0 ;
                0 0 0 0 0 v u 0 0 w 0 0 0 0 0 ;
                0 0 0 0 0 0 v u 0 0 w 0 0 0 0 ;
                0 0 0 0 0 0 0 v u 0 0 w 0 0 0 ;
                0 0 0 0 0 0 0 0 0 v u 0 w 0 0 ;
                0 0 0 0 0 0 0 0 0 0 v u 0 w 0 ;
                0 0 0 0 0 0 0 0 0 0 0 0 v u w ];

    limitpoint = decast * tobezier;
    limitpoint = decast(5:10, 6:15)  * limitpoint;
    
    if nargout > 2
        du             = 3 * [ limitpoint(2, :) - limitpoint(4, :) ;
                               limitpoint(3, :) - limitpoint(5, :) ;
                               limitpoint(5, :) - limitpoint(6, :) ];
        dv             = 3 * [ limitpoint(1, :) - limitpoint(4, :) ;
                               limitpoint(2, :) - limitpoint(5, :) ;
                               limitpoint(4, :) - limitpoint(6, :) ];
        
        dudu           = 4 * ( du(2, :) - du(3, :) );
        dudv           = 4 * ( du(1, :) - du(3, :) );
        dvdv           = 4 * ( dv(1, :) - dv(3, :) );
        
        secderiv       = [ dudu ; dudv ; dvdv ];
    end
    
    limitpoint = decast(8:10, 10:15) * limitpoint;
    deriv      = 4 * [ limitpoint(3, :) - limitpoint(1, :) ;
                       limitpoint(2, :) - limitpoint(3, :) ;
                       limitpoint(1, :) - limitpoint(2, :) ];
    limitpoint = decast(10, 13:15)   * limitpoint;

end

