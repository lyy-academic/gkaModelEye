function sceneGeometry = addEyePosePolyModel(sceneGeometry, varargin)
% Generate a polynomial model relating ellipse params to eyePose
%
% Syntax:
%  sceneGeometry = addEyePosePolyModel(sceneGeometry, varargin)
%
% Description:
%	This routine computes the pupil ellipse that is produce by a range of
%	eyePose values. The relationship between pupil ellipse values and
%	eyePose values are fit with a high-order polynomial model. This model
%	is then used to provide an initial x0 value for searches performed in
%	eyePoseEllipse fit and pupilProjection_inv.
%
% Inputs:
%   sceneGeometry         - Structure. SEE: createSceneGeometry
%
% Optional key/value pairs:
%  'eyePoseLB/UB'         - A 1x4 vector that provides the lower (upper)
%                           bounds on the eyePoses for the search [azimuth,
%                           elevation, torsion, stop radius].
%  'gridDensity'          - Scalar. How densely the eyePose grid of values
%                           is sampled. Computation time will grow
%                           geometrically with grid density.
%  'estimateCompTime'     - Logical. If set to true, the routine will exit
%                           without performing the computation and will
%                           instead report the expected time to compute the
%                           grid for these parameter settings.
%
% Outputs:
%   sceneGeometry         - Structure. The input structure is copied over
%                           to the output, with the addition of the field
%                           polyModel.
%
% Examples:
%{
    sceneGeometry = createSceneGeometry();
    sceneGeometry = calcEyePoseGrid(sceneGeometry);
%}
%{
    % Estimate computation time
    sceneGeometry = createSceneGeometry();
    calcEyePoseGrid(sceneGeometry,'gridDensity',85,'estimateCompTime',true);
%}

%% Parse input
p = inputParser;

% Required input
p.addRequired('sceneGeometry',@isstruct);

% Optional params
p.addParameter('eyePoseLB',[-89,-89,0,0.1],@isnumeric);
p.addParameter('eyePoseUB',[89,89,0,4],@isnumeric);
p.addParameter('gridDensity',25,@isscalar);
p.addParameter('polyModelOrder',6,@isscalar);
p.addParameter('estimateCompTime',false,@islogical);
p.addParameter('verbose', false, @islogical);

% Parse and check the parameters
p.parse(sceneGeometry, varargin{:});

% Anonymous function to produce exponentially spaced values
expVal = 20;
expspace = @(mn,mx,n) (mx-mn)/expVal * (10.^(linspace(0, log10(expVal+1), n)) - 1) + mn;

% Generate the grid of values
aziVals = unique([fliplr(-expspace(0,abs(p.Results.eyePoseLB(1)),round(p.Results.gridDensity/2))) expspace(0,p.Results.eyePoseUB(1),round(p.Results.gridDensity/2))]);
eleVals = unique([fliplr(-expspace(0,abs(p.Results.eyePoseLB(2)),round(p.Results.gridDensity/2))) expspace(0,p.Results.eyePoseUB(2),round(p.Results.gridDensity/2))]);
stopVals = linspace(p.Results.eyePoseLB(4),p.Results.eyePoseUB(4),round(p.Results.gridDensity/4));

nComps = length(aziVals)*length(eleVals)*length(stopVals);

% Check if we are just to estimate the computation time
if p.Results.estimateCompTime || p.Results.verbose
    tic
    for ii=1:5
        pupilProjection_fwd([0 0 0 1], sceneGeometry, 'nStopPerimPoints', 5);
    end
    compTimeMins = toc*nComps/5/60;
    fprintf('Grid computation time for these parameters is %2.2f mins\n',compTimeMins);
end

if p.Results.estimateCompTime
    fprintf('Returning without performing the computation\n');
    return
end

% Loop through the eyePoses and calculate the forward projection
eyePoses = [];
pupilEllipses = [];
tic
hold on
idx=1;
for xx=1:length(aziVals)
    % Update the progress display
    if p.Results.verbose && mod(idx,round(nComps/50))==0
        fprintf('\b.\n');
    end    
    for yy=1:length(eleVals)
        for zz=1:length(stopVals)
            eyePose = [aziVals(xx) eleVals(yy) 0 stopVals(zz)];
            pupilEllipse = pupilProjection_fwd(eyePose, sceneGeometry, 'nStopPerimPoints', 5);
            if ~any(isnan(pupilEllipse))
                eyePoses(idx,:) = eyePose;
                pupilEllipses(idx,:) = pupilEllipse;
                idx=idx+1;
            end
        end
    end
    toc
end
compTimeMins = toc/60;

% report completion of prediction grid
if p.Results.verbose
    fprintf('Finished grid construction. Now fitting polynomial model\n');
end

% Fit the polynomial model that relates pupilEllipse parameters to each of
% the eye pose parameters
sceneGeometry.polyModel.azimuth = ...
    polyfitn( sceneGeometry.eyePoseGrid.pupilEllipses, ...
    sceneGeometry.eyePoseGrid.eyePoses(:,1),p.Results.polyModelOrder);
sceneGeometry.polyModel.elevation = ...
    polyfitn( sceneGeometry.eyePoseGrid.pupilEllipses, ...
    sceneGeometry.eyePoseGrid.eyePoses(:,2),p.Results.polyModelOrder);
sceneGeometry.polyModel.stopRadius = ...
    polyfitn( sceneGeometry.eyePoseGrid.pupilEllipses, ...
    sceneGeometry.eyePoseGrid.eyePoses(:,4),p.Results.polyModelOrder);


% Store the meta data
sceneGeometry.eyePoseGrid.meta = p.Results;
sceneGeometry.eyePoseGrid.compTimeMins = compTimeMins;

end % function



