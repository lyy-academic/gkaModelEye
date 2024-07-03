function [eyePoints, pointLabels,rayPath_1st,rayPath_4th] = addGlint_lyy_parallel(eyePoints,pointLabels,sceneGeometry,p,eyePose,parallelRays)
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
% for gg = 1:nGlints
    
%     % Perform the computation using the passed function handle
%     [virtualImageRay, initialRay, intersectError,rayPath] = ...
%         glintRayFunc(glintSourceWorld(:, gg), eyePose, args{:});
    
%     % If we have a good trace, add the glint point and label
%     if intersectError < rayTraceErrorThreshold
%         eyePoints = [eyePoints; virtualImageRay(1,:)];
%         glintLabel = sprintf('glint_%02d',gg);
%         pointLabels = [pointLabels; {glintLabel}];
%         outputRay=virtualImageRay;
%     end

% end

nRays=size(parallelRays,2);
Rstruc = struct('azi',nan(3,3),'ele',nan(3,3),'tor',nan(3,3),'empty',true);


opticalSystemFixRL=sceneGeometry.refraction.cameraToMedium.opticalSystem;
opticalSystemRot_1st=sceneGeometry.refraction.glint.opticalSystem;
opticalSystemFixLR=sceneGeometry.refraction.mediumToCamera.opticalSystem;

% opticalSystem for lens_back
[os_4th,~,~,~]=assembleOpticalSystem(sceneGeometry.eye,'surfaceSetName','cameraToRetina','skipMagCalc',true);
os_4th(49,:)=os_4th(50,:); % remove the last row, retina, by covering with nans
os_4th(48,end)=-os_4th(47,end);

[os_4th_return,~,~,~]=assembleOpticalSystem(sceneGeometry.eye,'surfaceSetName','retinaToCamera','skipMagCalc',true);
os_4th_return=os_4th_return(sum(isnan(os_4th_return),2)~=size(os_4th_return,2),:);
os_4th(48+1:48+size(os_4th_return,1)-3,:)=os_4th_return(4:end,:);

opticalSystemRot_4th=os_4th;

rayPath_1st={};
rayPath_4th={};

for rr = 1:nRays
    initialRay=parallelRays{rr};
    [outputRay_1st,rayPath_1st{rr}] = threeSystemTrace(initialRay', eyePose, sceneGeometry.eye.rotationCenters, opticalSystemFixRL, opticalSystemRot_1st, opticalSystemFixLR,Rstruc);
    [outputRay_4th,rayPath_4th{rr}] = threeSystemTrace(initialRay', eyePose, sceneGeometry.eye.rotationCenters, opticalSystemFixRL, opticalSystemRot_4th, opticalSystemFixLR,Rstruc);
    eyePoints = [eyePoints; outputRay_1st(1,:)];
    eyePoints = [eyePoints; outputRay_4th(1,:)];
    glintLabel = sprintf('glint_1st_%02d',rr);
    pointLabels = [pointLabels; {glintLabel}];
    glintLabel = sprintf('glint_4th_%02d',rr);
    pointLabels = [pointLabels; {glintLabel}];


    rayPath_1st{rr}(end+1,:)=outputRay_1st(1,:)+outputRay_1st(2,:)*50;
    rayPath_4th{rr}(end+1,:)=outputRay_4th(1,:)+outputRay_4th(2,:)*50;

    
end



end % addGlint



%% local functions

function [outputRayEyeCoord,rayPath] = threeSystemTrace(inputRayEyeCoord, eyePose, rotationCenters, opticalSystemFixRL, opticalSystemRot, opticalSystemFixLR, Rstruc)

    % Ray trace through the fixed system (typically, camera to the medium
    % adjacent to the eye)
    outputRayEyeCoord = rayTraceQuadrics(inputRayEyeCoord, opticalSystemFixRL)';
    
    % outputRayStack=outputRayEyeCoord;
    
    % Counter-rotate the outputRayEyeCoord by the eye pose. This places the
    % outputRayEyeCoord so that the eyeCoordTarget is in a position equivalent
    % to if the eye had rotated.
    outputRayEyeCoord = rotateEyeRay(outputRayEyeCoord, eyePose, rotationCenters, 'inverse', Rstruc);
    
    % Ray trace through the eye system that is subject to rotation
    [outputRayEyeCoord,rayPath] = rayTraceQuadrics(outputRayEyeCoord', opticalSystemRot);
    outputRayEyeCoord=outputRayEyeCoord';
    rayPath=rayPath';
    
    % Return the ray to the original eyeWorld coordinate frame, which is
    % aligned with the optical axis of the eye
    outputRayEyeCoord = rotateEyeRay(outputRayEyeCoord, eyePose, rotationCenters, 'forward', Rstruc);
    
    % for i=1:(size(rayPath,1)/2)
    %     rayPath(2*i-1:2*i,:)=rotateEyeRay(rayPath(2*i-1:2*i,:),eyePose,rotationCenters,'forward',Rstruc);
    % end
    
    % outputRayStack=[outputRayStack;rayPath];
    % outputRayStack=rayPath;
    
    % Ray trace back through the fixed system that is not subject to rotation
    outputRayEyeCoord = rayTraceQuadrics(outputRayEyeCoord', opticalSystemFixLR)';
    
    % outputRayStack=[outputRayStack;outputRayEyeCoord];
    
    end