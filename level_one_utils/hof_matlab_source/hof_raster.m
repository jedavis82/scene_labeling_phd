%--------------------------------------------------------------------------
% hof_raster.m
% version 1.01
% 5/18/11
% Author: Andrew Buck
%
% This is a wrapper function which calls the appropriate MEX-file to
% calculate the histograms of forces for raster data. Inputs are checked
% for validity before calling the appropriate MEX-file.  The MEX-files
% interface with a C++ wrapper class written by Ozy Sjahputera for
% calculating the histograms of forces using the original code by Pascal
% Matsakis.  The histograms of forces concept is based on the following
% work:
%
% [1] P. Matsakis, Relations spatiales structurelles et interpretation
%     d'images, Ph.D. dissertation, IRIT, Universite Paul Sabatier,
%     Tolouse, France, 1998.
% [2] P. Matsakis, L. Wendling, "A New Way to Represent the Relative
%     Position of Areal Objects", PAMI, vol. 21, no. 7, pp. 634-643, 1999.
%
% Syntax:
%   H = hof_raster(imgA, imgB, r)
%   H = hof_raster(imgA, imgB, 'hybrid')
%   H = hof_raster(..., 'PropertyName', PropertyValue, ...)
%
% Description:
%   H = hof_raster(imgA, imgB, r) assigns H the F-r histogram between imgA
%   and imgB, where imgA and imgB are image matricies of the same dimension
%   in the range [0, 255] and r is the type of force (0 = constant, 2 =
%   gravitational).  imgA represents the argument object and imgB
%   represents the referent object.  A value of 0 means that the pixel does
%   not belong at all to the object, and the value 255 means that it
%   totally belongs to the object.  The object is fuzzy if there exists a
%   pixel whose gray level g is such that 0<g<255.  Otherwise it is crisp.
%
%   H = hof_raster(imgA, imgB, 'hybrid') assigns H the F-hybrid histogram
%   between imgA and imgB.
%
%   H = hof_raster(..., 'PropertyName', PropertyValue, ...) can be used
%   to set additional properties using either of the two previous methods. 
%   Valid properties are described below.
%
% Properties:
%   'NumberDirections' | (default = 180) : A positive multiple of 4 (e.g.,
%       16, 32, 64, 120, 180, 360).  Forces will be considered in
%       'NumberDirections' directions.
%
%   'FuzzyMethod' | (default = 1) : Method for computation of histograms
%       associated with fuzzy objects.  Values can be:
%            1 : simple sum scheme, see Krishnapuram, Keller and Ma, 1993
%            2 : double sum scheme, see Dubois and Jaulent, 1987
%           -1 : equivalent to the simple sum scheme, but pairs of pixels
%                are processed instead of alpha-cuts
%           -2 : equivalent to the double sum scheme
%
%   'p0' | (default = 0.01) : For the computation of hybrid forces only.
%       See [1].  p0 is related to the size of the smallest object.  If the
%       objects are "not too close" (according to p0), only gravitational
%       forces are exerted.  We should have 0<p0<p1
%
%   'p1' | (default = 3) : For the computation of hybrid forces only.
%       See [1].  p1 is related to the size of the intersection.  When the
%       objects intersect, or are "very close", constant forces appear.  We
%       should have 0<p0<p1.
%
%   Updates:
%       1.00 02/11/11 : Original function
%       1.01 05/18/11 : Added ismatrix compatability function
%--------------------------------------------------------------------------

function H = hof_raster(varargin)

%Define default parameters
numDir = 180;
fMethod = 1;
p0 = 0.01;
p1 = 3;

%Check that we have a valid number of inputs
if nargin < 3
    error('Too few inputs');
elseif nargin > 11
    error('Too many inputs');
elseif ~mod(nargin,2)
    error('Invalid number of parameters');
end

%Validate first three terms
imgA = varargin{1};
imgB = varargin{2};
htype = varargin{3};

%Check that images are valid
if ~ismatrix(imgA) || ~ismatrix(imgB)
    error('Images must be in matrix form');
elseif size(imgA,1) ~= size(imgB,1) || size(imgA,2) ~= size(imgB,2)
    error('Images must be the same size');
elseif ~isreal(imgA) || ~isreal(imgB) || ...
       max(max(~isfinite(imgA))) || max(max(~isfinite(imgB)))
    error('Images must contain real values in the range [0,255]');
elseif min(imgA(:)) < 0 || min(imgB(:)) < 0 || ...
       max(imgA(:)) > 255 || max(imgB(:)) > 255
    error('Images must contain real values in the range [0,255]');
elseif max(imgA(:)) == 0 || max(imgB(:)) == 0
    error('Images must contain at least one non-zero pixel');
end

%Determine if images are crisp or fuzzy
if max(max(imgA > 0 & imgA < 255)) || max(max(imgB > 0 & imgB < 255))
    fuzzy = 1;
else
    fuzzy = 0;
end  

%Convert images to uint8 format and transpose
%Transpose is required because MATLAB stores arrays by columns, while the
%C++ code requires the images to be stored by rows.
imgA = uint8(imgA');
imgB = uint8(imgB');

%Check for a valid 3rd term
if ischar(htype)
    if ~strcmpi(htype, 'hybrid')
        error('Histogram type must be a finite real number or ''hybrid''');
    end
elseif ~isscalar(htype)
    error('Histogram type must be a finite real number or ''hybrid''');
elseif ~(isreal(htype) && isfinite(htype))
    error('Histogram type must be a finite real number or ''hybrid''');
end

%Get parameter values
for i = 4:2:nargin
    paramName = varargin{i};
    if strcmpi(paramName, 'NumberDirections')
        numDir = varargin{i+1};
        if ~isscalar(numDir)
            error('''NumberDirections'' must be a positive finite real multiple of 4');
        elseif~(isreal(numDir) && isfinite(numDir))
            error('''NumberDirections'' must be a positive finite real multiple of 4');
        elseif numDir <= 0 || mod(numDir,4)
            error('''NumberDirections'' must be a positive finite real multiple of 4');
        end
    elseif strcmpi(paramName, 'FuzzyMethod')
        fMethod = varargin{i+1};
        if ~isscalar(fMethod)
            error('''FuzzyMethod'' must be either -2, -1, 1, or 2');
        elseif ~(isreal(fMethod) && isfinite(fMethod))
            error('''FuzzyMethod'' must be either -2, -1, 1, or 2');
        elseif ~max(fMethod == [-2 -1 1 2])
            error('''FuzzyMethod'' must be either -2, -1, 1, or 2');
        end
    elseif strcmpi(paramName, 'p0')
        p0 = varargin{i+1};
        if ~isscalar(p0)
            error('''p0'' must be a positive finite real value');
        elseif ~(isreal(p0) && isfinite(p0))
            error('''p0'' must be a positive finite real value');
        elseif p0 <= 0
            error('''p0'' must be a positive finite real value');
        end
    elseif strcmpi(paramName, 'p1')
        p1 = varargin{i+1};
        if ~isscalar(p1)
            error('''p1'' must be a positive finite real value');
        elseif ~(isreal(p1) && isfinite(p1))
            error('''p1'' must be a positive finite real value');
        elseif p1 <= 0
            error('''p1'' must be a positive finite real value');
        end
    end
end

%Make sure p0 and p1 are valid
if p0 >= p1
    error('''p1'' must be greater than ''p0''');
end

%Call the appropriate MEX-file
if strcmp(strtok(version), '7.9.0.529')
    if fuzzy
        if strcmpi(htype, 'hybrid')
            H = F02Histogram_FuzzyRaster_09(imgA, imgB, numDir, fMethod, p0, p1);
        else
            H = FRHistogram_FuzzyRaster_09(imgA, imgB, htype, numDir, fMethod);
        end
    else
        if strcmpi(htype, 'hybrid')
            H = F02Histogram_CrispRaster_09(imgA, imgB, numDir, p0, p1);
        else
            H = FRHistogram_CrispRaster_09(imgA, imgB, htype, numDir);
        end
    end
else
    if fuzzy
        if strcmpi(htype, 'hybrid')
            H = F02Histogram_FuzzyRaster(imgA, imgB, numDir, fMethod, p0, p1);
        else
            H = FRHistogram_FuzzyRaster(imgA, imgB, htype, numDir, fMethod);
        end
    else
        if strcmpi(htype, 'hybrid')
            H = F02Histogram_CrispRaster(imgA, imgB, numDir, p0, p1);
        else
            H = FRHistogram_CrispRaster(imgA, imgB, htype, numDir);
        end
    end
end


%Compatability function to determine if input is a matrix
function ism = ismatrix (u)

ism = (ndims(u) == 2) & (min(size(u)) ~= 1);