function coordinates = surfaceGrid(S, boundingBox, vertexDensity, polarGrid, bbTol)
% Returns a set coordinates evenly distributed on the quadric surface
%
% Syntax:
%  coordinates = qaudric.surfaceGrid(S, boundingBox, vertexDensity, polarGrid, bbTol)
%
% Description:
%   Returns a set of coordinates that are spaced across the quadric
%   surface. If a polarGrid is requested, then the points are evenly spaced
%   in geodetic coordinates. Otherwise, the points are evenly spaced in
%   Cartesian coordinates.
%
% Inputs:
%   S                     - 1x10 vector or 4x4 matrix of the quadric
%                           surface.
%   boundingBox           - 1x6 vector that specifies:
%                           	[xmin, xmax, ymin, ymax, zmin, zmax]
%                           These values set the bounds within which the
%                           coordinates are reported.
%   vertexDensity         - Scalar. The density of the mesh grid. Defaults
%                           to unity.
%   polarGrid             - Logical. If set to true, produces a polar
%                           distribution of points. This is valid only if
%                           the quadric is an ellipsoid.
%   bbTol                 - Scalar. Defines the tolerance within which the
%                           intersection must be within the boundingBox.
%                           Default value is 0.1. Handles the situation in
%                           which the intersection is right on the boundary
%                           but is numerically outside.
%
% Outputs:
%   coordinates           - nx3 matrix that provides the [x; y; z]
%                           coordinates of n points that are on the surface
%                           of the quadric within the boundingBox
%
% Examples:
%{
    S = quadric.scale(quadric.unitSphere,[4 5 3]);
    boundingBox = [0 50 -30 30 -20 20];
    coordinates = quadric.surfaceGrid(S,boundingBox,50);
    plot3(coordinates(:,1),coordinates(:,2),coordinates(:,3),'.r')
%}

% Handle incomplete input arguments
if nargin == 2
    vertexDensity = 1;
end

if nargin == 3
    polarGrid = true;
    bbTol = 1e-2;
end

if nargin==4
    bbTol = 1e-2;
end

% Obtain the surfaceGrid
if polarGrid
    % Assemble the set of coordinates
    coordinates = [];
    for latitude = linspace(-180,180,vertexDensity*2)
        for longitude=linspace(-180,180,vertexDensity)
            X = quadric.geodeticToCart( [latitude; longitude; 0], S );
            % Store the coordinate value.
            coordinates(end+1,:)= X;
        end
    end
    % Remove the coordinates that are outside the bounding box
    retainCoords = (coordinates(:,1) > boundingBox(1)-bbTol) .* (coordinates(:,1) < boundingBox(2)+bbTol) .* ...
        (coordinates(:,2) > boundingBox(3)-bbTol) .* (coordinates(:,2) < boundingBox(4)+bbTol) .* ...
        (coordinates(:,3) > boundingBox(5)-bbTol) .* (coordinates(:,3) < boundingBox(6)+bbTol);
    coordinates = coordinates(logical(retainCoords),:);
else
    % Produce a set of coordinates that are linearly spaced across the
    % boundingBox.
    coordinates = [linspace(boundingBox(1),boundingBox(2),vertexDensity); ...
	linspace(boundingBox(3),boundingBox(4),vertexDensity); ...
    	linspace(boundingBox(5),boundingBox(6),vertexDensity)];
end

end % surfaceGrid

