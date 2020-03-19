function [headPoints, pointLabels] = applyEyeRotation(eyePoints,pointLabels,sceneGeometry,p,eyePose)
% Subject the eyeWorld points to the effect of eye rotation
%
% Syntax:
%  [headPoints, pointLabels] = applyEyeRotation(eyePoints,pointLabels,sceneGeometry,p,eyePose)
%
% Description:
%   Rotate the eyeWorld points so that they are in the position specified
%   by the eyePose. As the eyeWorld coordinate frame is defined with
%   respect to the optical axis of the eye, the matrix of points that
%   result from this operation are termed as being in "head world"
%   coordinates, which is the coordinate frame defined with respect to an
%   eye that is locked in the [0 0 0] eye pose.
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
%   headPoints            - nx3 vector. Points in head world coordinates.
%   pointLabels           - nx1 cell array. The name of each eye point.
%


% Omit the eye rotation centers themselves from rotation.
rotatePointsIdx = find(~contains(pointLabels,{'Rotation'}));

% Copy the eyePoints to the headPoints
headPoints = eyePoints;

% Loop through the points to be rotated. We pass the rotation matrix to
% avoid having to re-calculate this for the rotation of each point.
R = [];
for pp = 1:length(rotatePointsIdx)
    headPoints(rotatePointsIdx(pp),:) = rotateEyeCoord(...
        eyePoints(rotatePointsIdx(pp),:), ...
        eyePose, ...
        sceneGeometry.eye.rotationCenters, ...
        'forward', ...
        R);
end

% If we are projecting a full eye model, label as hidden those posterior
% segment points that are posterior to the most posterior of the centers of
% rotation of the eye, and thus would not be visible to the camera.
if p.Results.fullEyeModelFlag
    seenIdx = strcmp(pointLabels,'retina') .* (headPoints(:,1) >= min([sceneGeometry.eye.rotationCenters.azi(1) sceneGeometry.eye.rotationCenters.ele(1)]));
    seenIdx = logical(seenIdx + ~strcmp(pointLabels,'retina'));
    pointLabels(~seenIdx) = strcat(pointLabels(~seenIdx),'_hidden');
end

end % applyEyeRotation