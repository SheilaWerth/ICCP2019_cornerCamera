function [A,removeRows] = removeCorner(A,numRemove)
% this function removes a square pixel swath (length numRemove) from the 
% floor corner
% removeRows is the indices of the rows removed
imgDim = length(A(1,:));
removeMask = zeros(imgDim,imgDim);
removeMask(1:numRemove,:) = 1;

indexMask = reshape(1:imgDim^2,[imgDim,imgDim]);

removeMask = removeMask.*indexMask;
removeRows = find(removeMask); 

A(removeRows,:) = [];

% and remove columns too - we can't reconstruct these angles without the
% data that we removed
% A(:,end-numRemove+1:end) = [];


end

