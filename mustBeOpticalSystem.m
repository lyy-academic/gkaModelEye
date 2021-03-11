function mustBeOpticalSystem(a)
% Argument validation function that requires optical system, eye, or scene
%
% Syntax:
%  mustBeOpticalSystem(a)
%
% Description:
%   For use in the arguments validation block in the calling function
%

validFlag = false;

% If we were passed a structure, ensure that it has the components of
% either a sceneGeometry, eye structure, or system set
if isstruct(a)
    if isfield(a,'cornea') || isfield(a,'eye') || isfield(a,'opticalSystem')
        validFlag = true;
    end
end

% If this is a matrix, then check that it has 19 elements in the second
% dimension
if ismatrix(a)
    if size(a,2) == 19
        validFlag = true;
    end
end

% Throw an error if not valid
if ~validFlag
    eidType = 'mustBeOpticalSystem:notValidMatrixOrStruct';
    msgType = 'Input must be either an optical system matrix, or an eye or sceneGeometry structure.';
    throwAsCaller(MException(eidType,msgType))
end