function [retinaPoints,rayPaths] = calcRayBundleFromField(opticalSystem,fieldAngularPosition,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord,rayIntersectionHeight,cameraMedium)
% Returns a nodal ray that arises from the specified field location
%
% Syntax:
%  [rayPath,angleError] = calcNodalRayFromField(opticalSystem,fieldAngularPosition,rayOriginDistance,angleReferenceCoord,cameraMedium)
%
% Description
%   
%
% Inputs:
%   opticalSystem         - Either an eye structure (from which a
%                           "mediumToRetina" optical system in air will be
%                           derived), or an mx19 matrix, where m is set by
%                           the key value opticalSystemNumRows. Each row
%                           contains the values:
%                               [S side bb must n]
%                           where:
%                               S     - 1x10 quadric surface vector
%                               side  - Scalar taking the value -1 or 1
%                                       that defines which of the two
%                                       points of intersection on the
%                                       quadric should be used as the
%                                       refractive surface.
%                               bb    - 1x6 vector defining the bounding
%                                       box within which the refractive
%                                       surface is present.
%                               must  - Scalar taking the value of 0 or 1,
%                                       where 1 indicates that the ray must
%                                       intersect the surface. If the ray
%                                       misses a required surface, the
%                                       routine exits with nans for the
%                                       outputRay.
%                               n     - Refractive index of the surface.
%   fieldAngularPosition  - 2x1 vector that provides the coordinates of the
%                           origin of the nodal ray in [horizontal,
%                           vertical[ degrees with respect to the
%                           coordinate specified in angleReferenceCoord.
%   rayOriginDistance     - Scalar. The distance (in mm) of the origin of
%                           the ray from the distanceReferenceCoord.
%   angleReferenceCoord   - 3x1 vector that provides the coordinate from
%                           which the ray origin angles and distance are
%                           to be calculated. By default, this is [0;0;0],
%                           which is the origin coordinate on the
%                           longitudinal axis.
%   distanceReferenceCoord - 3x1 vector that provides the coordinate from
%                           which the rayOriginDistance is calculated. The
%                           The principal point is a typical choice. If not
%                           defined, is set to [0;0;0], which is the origin
%                           coordinate on the longitudinal axis.
%   cameraMedium          - String. The medium in which the eye is located.
%                           Defaults to 'air'.
%
% Outputs:
%   rayPath               - 3xm matrix that provides the ray coordinates
%                           at each surface. The value for rayPath(1,:)
%                           is equal to initial position. If a surface is
%                           missed, then the coordinates for that surface
%                           will be nan.
%   angleError            - Scalar. The departure from parallel of the 
%                           incident and emergent rays (deg)
%
% Examples:
%{
    % Define a default model eye
    accommodation = 1;
    eye = modelEyeParameters('accommodation',accommodation);
    % Pick a field location
    rayOriginDistance = 1000/accommodation;
    angleReferenceCoord = eye.landmarks.incidentNode.coords;
    distanceReferenceCoord = calcPrincipalPoint(eye);
    fieldAngularPosition = [10 10];
    % Trace the bundle
    retinaPoints = calcRayBundleFromField(eye,fieldAngularPosition,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord);
    % Confirm that the angleError is within tolerance
    
%}


arguments
    opticalSystem
    fieldAngularPosition (1,2) {mustBeNumeric} = [0,0]
    rayOriginDistance (1,1) {mustBeNumeric} = 1500
    angleReferenceCoord (3,1) {mustBeNumeric} = [0;0;0]
    distanceReferenceCoord (3,1) {mustBeNumeric} = [0;0;0]
    rayIntersectionHeight (1,1) {mustBeNumeric} = 1
    cameraMedium = 'air'
end


% Check if we were passed an eye model. If so, create the optical system
if isstruct(opticalSystem)
    if isfield(opticalSystem,'cornea')
        eye = opticalSystem;
        clear opticalSystem;
        opticalSystem = assembleOpticalSystem(eye,...
            'surfaceSetName','mediumToRetina','cameraMedium',cameraMedium);
    end
end


% Define the primary fieldRay
fieldRay = calcFieldRay(fieldAngularPosition,rayOriginDistance,angleReferenceCoord,distanceReferenceCoord);

% Reverse the ray direction
fieldRay(:,2) = -fieldRay(:,2);

% Decompose the fieldRay
rayOrigin = fieldRay(:,1);
[p1p2, p1p3] = quadric.rayToAngles(fieldRay);

% Loop over a set of displacement angles
retinaPoints = [];
rayPaths = {};
bundleAngle = atand(rayIntersectionHeight/rayOriginDistance);
nDivisions = 5;
for hh = linspace(-bundleAngle,bundleAngle,nDivisions)
    for vv = linspace(-bundleAngle,bundleAngle,nDivisions)
        if norm([hh vv])>bundleAngle
            continue
        end
        bundleRay = quadric.anglesToRay(rayOrigin,p1p2+hh,p1p3+vv);
        [outputRay, rayPaths{end+1}] = rayTraceQuadrics(bundleRay, opticalSystem);
        retinaPoints(:,end+1) = outputRay(:,1);
    end
end


end