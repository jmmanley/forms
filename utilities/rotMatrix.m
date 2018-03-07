function R = rotMatrix(x,y,z)
% Jason Manley, Sep 2017

Rx = [1 0 0 0;
    0 cos(x) -sin(x) 0;
    0 sin(x) cos(x) 0;
    0 0 0 1];

Ry = [cos(y) 0 sin(y) 0;
    0 1 0 0;
    -sin(y) 0 cos(y) 0;
    0 0 0 1];

Rz = [cos(z) -sin(z) 0 0;
    sin(z) cos(z) 0 0;
    0 0 1 0;
    0 0 0 1];

R = Rx*Ry*Rz;


end