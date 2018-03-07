function [ points3d points2d onsil ] = parse_fpjconstraints(string)

% PARSE_FPJCONSTRAINTS  Parse the syntax used for point constraints
%   in .FPJ files
%
%   See also READ_FPJ, READ_PROJECT

if isempty(string)
    points3d = zeros(0, 1);
    points2d = zeros(0, 2);
    onsil    = zeros(0, 1);
else
    C = textscan(string, '%u:%f:%f:%u', 'delimiter', ',');

    points3d = C{1} + 1;
    points2d = [C{2} C{3}];
    onsil    = C{4};
end

end

