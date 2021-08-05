function [ freqWeight, onedbasisfns, A_penumbra, removeRows] = ...
    formCmat_dc_cornerCrop( m,n,w,S,pow, numRemove, dctFlag, dctType)
%FORMATAMAT this function constucts the c_synthesis matrix, [B A] that
%operates on the Z vector [x; u], where x are the coefficients of the
%penumbra and u are the coefficients of the floor
% m is the number of pixels in the image
% n is the number of coefficients in 
% freqWeight is the 1/f prior. Diagonal matrix
% basisFns - inverse fourier transform
%
% first make A_synthesis! This is the matrix that operates on teh floor
% Generate identity matrix
identity = eye(n);
% initialize the format of output matrix M
A = zeros(m,n);

% % loop through columns of identity
for Icol = 1:n
    % input into transform
    thisCol = identity(:,Icol);
    % perform decomposition
    A(:,Icol) = A_synthesis( thisCol,w,S);   
    
end

% and now generate B matrix - that operates on the penumbra
onedbasisfns = [];

%%
if dctFlag == 0 % use Fourier basis
    
    for yfreq = 1:sqrt(m)-numRemove
        onedfreqim = zeros(sqrt(m),1);
        onedfreqim(yfreq) = 1;
        onedim = ifft(ifftshift(onedfreqim));
        tmp5 = real(onedim);
        tmp6 = imag(onedim);
        onedbasisfns = [onedbasisfns, tmp5(:), tmp6(:)];
    end
    onedbasisfns = onedbasisfns(:, [1:(sqrt(m))]);
    onedbasisfns2 = onedbasisfns(:,1:end-1)./norm(onedbasisfns(:,1),2);
    % onedbasisfns2(:,31) = onedbasisfns(:,end)./norm(onedbasisfns(:,end),2);
    % size(onedbasisfns)
    onedbasisfns2(:,end+1) = onedbasisfns(:,end)./norm(onedbasisfns(:,end),2);
    onedbasisfns = onedbasisfns2;
    
else % use the DCT
    for yfreq = 1:sqrt(m)
        onedfreqim = zeros(sqrt(m),1);
        onedfreqim(yfreq) = 1;
        onedim = idct(onedfreqim,'Type', dctType);
        onedim = onedim/norm(onedim,2);
        onedbasisfns = [onedbasisfns, onedim];
    end
    onedbasisfns = fliplr(onedbasisfns);
end
%%
% sine and cosine frequencies
tmp = (sqrt(m)-1)/2:-1:1 ;
tmp2 = repmat(tmp,[2,1]);
onedfreqs = tmp2(:);
onedfreqs = [onedfreqs; 1];

% but we still need to convert this to an image rendered on the ground.
freqWeight = diag(1./onedfreqs.^pow);
tmpB = onedbasisfns * freqWeight;

params = cornerParams(sqrt(m)); % input the img dim
A_penumbra= getAmat(params, ones([1,sqrt(m),3]));
A_penumbra = reshape(A_penumbra,[sqrt(m)*sqrt(m),sqrt(m)]);
[A_penumbra,removeRows] = removeCorner(A_penumbra,numRemove);




end

