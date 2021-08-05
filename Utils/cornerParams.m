function [ params ] = cornerParams( M )
% initializes corner rendering parameters
%   For the static corner camera readout version.
% modified from corner/initParams.m, % dec. 8, 2017 billf created.
% July 3, 2017  billf created.


% inSceneH, V is the horiz and vert pixel count of the input scene.
% estSceneH, V is same, but for the estimated scene (we'll estimate
% just a 1-d scene, for computational efficiency).
% imageH, V  is the horiz and vert pixel count for the rendered image
% on the ground near the corner.  See explanation.key for a figure.
% imageCorner1,2,3,4 are the x,y,z positions of each of the four corner
% points of the image plane we record from.
% imageNormal is a 3-vector giving a surface normal of the image plane.
% vectorToObserver is a 3-vector pointing to the observer/camera.

% a 2-d image projecting onto the corner
%%%% 'amatBlurKernel', [1 4 6 4 1]/16, ...
params = struct(...
    'sceneToCorner', 20, ...
    'sceneVerticalOffset', 10, ...
    'sceneRenderThetaMax', pi/2, ...
    'sceneRenderDeltaTheta', pi/2, ...
    'imageH', M, ...
    'imageV', M, ...
    'imageCorner1', [-1;0;0], ...
    'imageCorner2', [0;0;0], ...
    'imageCorner3', [0; -1; 0], ...
    'imageCorner4', [-1; -1; 0], ...
    'imageNormal', [0;0;1], ...
    'brdf', 'Lambertian', ...
    'amatBlurKernel', [1 8 28 56 70 56 28 8 1]/256, ...
    'amatUpsampleFactor', 2, ...
    'unitVectorToObserver', [0;0;1]); % assuming we are looking straight down at the floor


