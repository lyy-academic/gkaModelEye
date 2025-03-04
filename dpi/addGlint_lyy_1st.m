function [eyePoints, pointLabels,initialRay, rayPath,outputRay] = addGlint_lyy_1st(eyePoints,pointLabels,sceneGeometry,p,eyePose)
% Add glint(s) to the eye model
%
% Syntax:
%  [eyePoints, pointLabels] = addGlint(eyePoints,pointLabels,sceneGeometry,p,eyePose)
%
% Description:
%	Find the location of the glint in the eye world coordinate frame. The
%	glint is the reflection of a light source from the tear film of the
%	eye. The location of the glint in the image is subject to refraction by
%   artificial lenses.
%
% Inputs:
%   eyePoints             - nx3 vector. Points in eye world coordinates.
%   pointLabels           - nx1 cell array. The name of each eye point.
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%   p                     - Structure. The structure returned by the
%                           parameter parser in the calling function.
%   eyePose               - A 1x4 vector provides values for [eyeAzimuth,
%                           eyeElevation, eyeTorsion, stopRadius].
%                           Azimuth, elevation, and torsion are in units of
%                           head-centered (extrinsic) degrees, and stop
%                           radius in mm.
%
% Outputs:
%   eyePoints             - nx3 vector. Points in eye world coordinates.
%   pointLabels           - nx1 cell array. The name of each eye point.
%


% Extract some values for clarity in the code that follows
glintRayFunc = p.Results.glintRayFunc;
rayTraceErrorThreshold = p.Results.rayTraceErrorThreshold;

% If we do not have the sceneGeometry components needed for the
% calculation, return.
if ~isfield(sceneGeometry,'refraction') || isempty(glintRayFunc)
    return
end
if ~isfield(sceneGeometry.refraction,'glint')
    return
end

% Get the position of the light source(s) in the world coordinate frame
% that is the source of the glint
glintSourceWorld = sceneGeometry.cameraPosition.translation + ...
    sceneGeometry.cameraPosition.glintSourceRelative;

% How many light sources do we have?
nGlints = size(glintSourceWorld,2);

% Assemble the args
args = {sceneGeometry.cameraPosition.translation, ...
    sceneGeometry.eye.rotationCenters, ...
    sceneGeometry.refraction.cameraToMedium.opticalSystem, ...
    sceneGeometry.refraction.glint.opticalSystem, ...
    sceneGeometry.refraction.mediumToCamera.opticalSystem};

% Loop through the glints
for gg = 1:nGlints
    
    % Perform the computation using the passed function handle
    [virtualImageRay, initialRay, intersectError,rayPath] = ...
        glintRayFunc(glintSourceWorld(:, gg), eyePose, args{:});
    
    % If we have a good trace, add the glint point and label
    if intersectError < rayTraceErrorThreshold
        eyePoints = [eyePoints; virtualImageRay(1,:)];
        glintLabel = sprintf('glint_%02d',gg);
        pointLabels = [pointLabels; {glintLabel}];
        outputRay=virtualImageRay;
    end
    
end


end % addGlint