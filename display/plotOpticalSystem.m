function figHandle = plotOpticalSystem(varargin)
% Creates a cross-section schematic illustration of the model eye
%
% Syntax:
%  plotModelEyeSchematic(eye)
%
% Description:
%   Create a schematic diagram of the model eye specified in the passed
%   eye variable.
%
% Inputs:
%   eye                   - An eye struct returned from
%                           modelEyeParameters()
%
% Optional key/value pairs:
%  'view'                 - String. The view to display. Valid choices
%                           include {'axial','sagittal'};
%  'newFigure'            - Logical. Determines if we create a new figure.
%  'plotColor'            - String. Matlab line spec code for line color,
%                           e.g., {'k','r','b'};
%
% Outputs:
%   figHandle             - Handle to a created figure. Empty if a new
%                           figure was not requested.
%
% Examples:
%{
    % Basic call for an axial and sagittal view plot in new windows
    eye = modelEyeParameters;
    plotModelEyeSchematic(eye);
    plotModelEyeSchematic(eye,'view','sag');
%}
%{
%}

%% input parser
p = inputParser; p.KeepUnmatched = true;

% Optional
p.addParameter('newFigure',true,@islogical);
p.addParameter('visible',true,@islogical);
p.addParameter('opticalSystem',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('rayPath',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('outputRay',[], @(x)(isempty(x) | isnumeric(x)));
p.addParameter('addLighting',false, @islogical);
p.addParameter('surfaceColors', ...
    {'',[0.5 0.5 0.5],[0.5 0.6 0.5],[0.5 0.7 0.5],[0.5 0.8 0.5],[0.5 0.9 0.5],[0.5 1.0 0.5],[0.5 1.0 0.5],[0.5 0.9 0.5],[0.5 0.8 0.5],[0.5 0.7 0.5],[0.5 0.6 0.5],[0.5 0.5 0.5],'blue','blue'},...
    @iscell);
p.addParameter('surfaceAlpha', ones(1,15)*0.1,@isnumeric);

p.addParameter('rayColor','red',@(x)(ischar(x) | isnumeric(x)));

% parse
p.parse(varargin{:})


% Open a figure
if p.Results.newFigure
    if p.Results.visible
        figHandle = figure('Visible', 'on');
    else
        figHandle = figure('Visible', 'off');
    end
else
    figHandle = gcf;
end
hold on

% Plot the opticalSystem if provided
if ~isempty(p.Results.opticalSystem)
    for ii=2:size(p.Results.opticalSystem,1)
        % Obtain a function handle for the polynomial
        F = quadric.vecToFunc(p.Results.opticalSystem(ii,1:10));
        boundingBox = p.Results.opticalSystem(ii,12:17);
        % Plot the surface
        plotSurface(F,boundingBox,p.Results.surfaceColors{ii},p.Results.surfaceAlpha(ii))
        hold on
    end
end

% Plot the rayPath if provided
if ~isempty(p.Results.rayPath)
    plot3(p.Results.rayPath(1,:),p.Results.rayPath(2,:),p.Results.rayPath(3,:),'-r');
end

% Add the outputRay if provided
if ~isempty(p.Results.outputRay)
    outputRay = p.Results.outputRay;
    p1=outputRay(:,1);
    p2=p1+outputRay(:,2).*3;
    r = [p1 p2];
    plot3(r(1,:),r(2,:),r(3,:),'-g');
end

% Add a lighting source if requested
if p.Results.addLighting
    camlight
    lighting gouraud
end

end % plotOpticalSystem


function plotSurface(F,boundingBox,surfColor,surfAlpha)

[xx, yy, zz]=meshgrid( linspace(boundingBox(1),boundingBox(2),100),...
    linspace(boundingBox(3),boundingBox(4),100),...
    linspace(boundingBox(5),boundingBox(6),100));
    vertices = isosurface(xx, yy, zz, F(xx, yy, zz), 0);

p = patch(vertices);
p.FaceColor = surfColor;
p.EdgeColor = 'none';
alpha(surfAlpha);
daspect([1 1 1])
view(3); 
axis tight
axis equal

end