function [freqWeight, onedbasisfns] = newPriorWandD(m, dctType, pow)
onedbasisfns = [];

for yfreq = 1:m
    onedfreqim = zeros(m,1);
    onedfreqim(yfreq) = 1;
    onedim = idct(onedfreqim,'Type', dctType);
    onedim = onedim/norm(onedim,2);
    onedbasisfns = [onedbasisfns, onedim];
end
onedbasisfns = fliplr(onedbasisfns);

% sine and cosine frequencies
tmp = (m)/2:-1:1 ;
% size(tmp)
tmp2 = repmat(tmp,[2,1]);
% size(tmp2)
onedfreqs = tmp2(:);
% size(onedfreqs)
% onedfreqs = [onedfreqs; 1];

% but we still need to convert this to an image rendered on the ground.
freqWeight = diag(1./onedfreqs.^pow);

end

