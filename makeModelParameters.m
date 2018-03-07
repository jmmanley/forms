function parameters = makeModelParameters(parameters)
% Initializes the parameters for fitting a 3D morphable model to 2D images,
% as described in Cashman & Fitzgibbon, 2012, "What Shape Are Dolphins? 
% Building 3D Morphable Models from 2D Images"
%
% The default parameters are those utilized by Cashman & Fitzgibbon to fit
% the model to dolphins.
%
% J. Manley, updated 9/2017

if nargin < 1
    parameters = struct();
end

if ~isfield(parameters,'M')
    parameters.M = 8;
    %is the number of deformation basis shapes (D in the paper)
end

if ~isfield(parameters,'xi_0')
    parameters.xi_0 = 0.5;
    %is the weight on E^tp_0
end

if ~isfield(parameters,'xi_def')
    parameters.xi_def = 0.25;
    %is the weight on E^tp_m for m > 0
end

if ~isfield(parameters,'sigma_norm')
    parameters.sigma_norm = 10;
    %is the noise estimate on normals
end

if ~isfield(parameters,'sigma_con')
    parameters.sigma_con = 0.2;
    %is the noise estimate on point constraints
end

if ~isfield(parameters,'sigma_sil')
    parameters.sigma_sil = 1;
    %is the noise estimate on silhouettes
end


end