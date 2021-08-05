clear; 
clc; clear;close all

warning('off','all')

rmpath(genpath('/Users/sheilawerth/Dropbox/Lab/'))
addpath(genpath('/Users/sheilawerth/Dropbox/Lab/Passive corner camera code review/ICCP/Figure 5'))
rmpath(genpath('/Users/sheilawerth/Dropbox/Lab/Passive corner camera code review/ICCP/Figure 5/Utils/deleteMe'))




display = true;

%% PARAMETERS

pow = 1; % for frequency weighting
dctFlag = 1; % use dct or wavelet basis?
dctType = 2;
dctTypeNewPrior = 1;
freqWeightPowNewPrior = 1;
numRemove = 0; % how many rows to remove from the top of measurement to remove wall contribution
%note! If you remove n rows, disconsider the n greatest angles in our scene
%estimate.





%% LOAD DATA
datStr = 'figure7_blueChannel.mat';
load(datStr)


y1 = data.y1;
y2 = data.y2;
y = data.y;

M = length(y);

%% A MATRICES
imgDim = 93; % number of pixels along one edge of photograph, also the number of pixels in scene reconstruction

% decide how many levels of our wavelet transform we can compute
l = wmaxlev(imgDim,'db1');

% Wavelet transform
[coefs,coef_ordering] = wavedec2(zeros(imgDim,imgDim),l,'db1');

[f_weight, onedbasisfns, A_penumbra, removeRows ] =...
    formCmat_dc_cornerCrop( imgDim^2,length(coefs),'db1',...
    coef_ordering,pow, numRemove, dctFlag, dctType);
d = diag(f_weight);

% new prior gaussian part
[f_weight2,onedbasisfns2] = newPriorWandD(imgDim-1,dctTypeNewPrior,freqWeightPowNewPrior);
d2 = diag(f_weight2);


y(removeRows) = [];
y1(removeRows) = [];
y2(removeRows) = [];


penumbraIncident = data.penumbra_light; % just incident light, floor albedo is uniform here
penumbraIncident(removeRows) = [];
trueAmbientLight = data.ambient_light;
trueAmbientLight = median(trueAmbientLight);

if dctFlag == 1
    W =@(x) (onedbasisfns'*x)./d;
    W_adj =@(x) onedbasisfns*(x./d);
else % otherwise use wavelet basis
    a = ones(imgDim,1);
    
    level = fix(log2(imgDim));
    wavename = 'db4';
    [Coef,L] = wavedec(a,level,wavename);
    wavebasis = zeros(length(Coef),imgDim);
    for ii = 1:imgDim
        x = zeros(imgDim,1);
        x(ii) = 1;
        wavebasis(:,ii) = wavedec(x,level,wavename);
    end
    
    W =@(x) wavebasis*x;
    W_adj =@(theta) wavebasis'*theta;
end


%% MEASURED DATA and scaling
scaling = trueAmbientLight;
y1 = y1/scaling;
y = y/scaling;
f_truth = reshape(y1,[imgDim-numRemove,imgDim]);

S = 40 * A_penumbra; % A_fan
S = S/scaling;
ST = S';

%% Use white floor measurement to get cleaner reference, since we don't have the ground truth
lam2 = 10;

vp = zeros(imgDim,1);
u = vp;
tp = 1;
s = svd(S);
stepsize = 1/(s(1)^2 + lam2*max(1./(d.^2)));
for iter = 1:10000
    v = u - stepsize*(ST*(S*u - penumbraIncident/(trueAmbientLight)) + lam2*W_adj(W(u)));
     v(v<0) = 0;
    t = (sqrt(1+4*tp^2)+1)/2;
    u =  v + (tp-1)/t*(v-vp);
    if norm(v-vp)<1e-14
        break
    end
    tp = t;
    vp = v;   
end
sceneWhiteFloor = v;


figure;plot((1:imgDim)/imgDim*90,sceneWhiteFloor,'.-');
xlim([6,90-5])
drawnow 
penumbraWhiteFloor = S*v;
grid on
xlabel('angle')
ylabel('intensity')
title('scene estimate on white floor')

figure;
subplot(1,3,1)
imagesc(reshape(penumbraIncident/(trueAmbientLight)...
    ,[imgDim-numRemove,imgDim]));
axis image;colorbar;title('measured')
subplot(1,3,2)
imagesc(reshape(penumbraWhiteFloor...
    ,[imgDim-numRemove,imgDim]));
axis image;colorbar;title('estimte (white floor)')
subplot(1,3,3)
imagesc(reshape(penumbraWhiteFloor - penumbraIncident/(trueAmbientLight)...
    ,[imgDim-numRemove,imgDim]));
axis image;colorbar;title('error')
drawnow


%% Define Operators
N1 = length(y);

A = S; 


% For the iterative method
numOuterIter = 100; % number of times to alternate

dataCost = zeros(numOuterIter,1);


% tuning parameters
% 1. hidden scene:
lam2 = 1; % tuning parameter to promote smoothness of scene estimate

% 2. floor albedo:
lam5 = 1; % gaussian part of floor prior
lam4 = 0;%10 % floor prior tuning parameter - often times 0 is fine here. Turn up if needed 

%% sheila new prior stuff
% remove edges of the scene
P_para = S/(S'*S)*S';


numRemoveBothSidesScene = 1;

S3 =@(x) x(numRemoveBothSidesScene+1:end-numRemoveBothSidesScene);
S3_adj =@(x) [zeros(numRemoveBothSidesScene,1); x; zeros(numRemoveBothSidesScene,1)]; 

S4 = eye(imgDim);
S4 = S4(2:end,:);
S4 = sparse(S4);
S4_adj = S4';


q =S'*S + lam5*S4_adj*(onedbasisfns2*diag(d2).^2*onedbasisfns2'*S4);
qinv = inv(q);
storeMe = (q\(S'*P_para))';
storeMe2 = qinv*S'*P_para;
ZtZ =@(x)  storeMe * S3_adj(S3( storeMe2*x )) ;


storeMe3 = qinv*S'* P_para;
% 
% 
% 
paraCost =@(x) norm(S3( storeMe3 * x ),2)^2/2;
plotMe =@(x) ( storeMe3 * x ); 

plotMeTruth = plotMe([y1]);
plotMeAmbient = plotMe(data.ambient_light);
% projCostTruth = paraCost([onedbasisfns*(d.*x2_truth);y1]);

%%

singA = svd(A,'econ');
L = max(singA)^2 + lam2*max(1./d)^2;
stepsize_v = 1/L;

f_est = zeros(N1,1);
v_est = zeros(imgDim,1);%sceneWhiteFloor;

if dctFlag ==0
    p = zeros(length(Coef),1);
end


%%
maxIterSceneEstimate = 100;
maxIterFloorEstimate = 500;
stepsize_f = 0.0013; % floor stepsize

for outerIter = 1: numOuterIter
    
    
    % Estimate scene

    C =@(x) f_est.*(A*x);
    CT =@(x) A'*(f_est.*x);
    vp  = v_est;
    u = vp;
    tp = 1;
    
      if outerIter>1
        for iter = 1:maxIterSceneEstimate
            if dctFlag == 1
                v = u - stepsize_v*( CT(C(u) - (y - f_est)) + lam2*W_adj(W(u)) );
                v(v<0) = 0;
            else
                v = u - stepsize_v*(CT(C(u) - (y - f_est)));
                
                [v,p] = constrained_ell1(v, lam2, W, W_adj, 100, [0,inf],p);
            end
            t = (sqrt(1+4*tp^2)+1)/2;
            u =  v + (tp-1)/t*(v-vp);
            if norm(v-vp)<1e-8
                break
            end
            tp = t;
            vp = v;
        end
        v_est = v;
        coef_est = (onedbasisfns'*v_est);
      end
      
      
    % Estimate floor
    f_meas = y;
    B = (A*v_est+ones(N1,1));

    f = f_est;
    fp = f;
    st = f;
    t = 1;
    
    totalCost = zeros(1,maxIterFloorEstimate);
    dataCost = zeros(1,maxIterFloorEstimate);
    projCost = zeros(1,maxIterFloorEstimate);
    for iter = 1: maxIterFloorEstimate

        f = st - stepsize_f*(B.*(B.*f-f_meas) + lam4* storeMe * S3_adj(S3( storeMe2*f )) );
        
        
        f = f(:);
        f(f<0.05) = 0.05; f(f>1.0) = 1.0;
        
        dataCost(iter) = norm(f_meas - B.*f )^2/2;
        projCost(iter) = lam4*paraCost(f);
        totalCost(iter) = dataCost(iter) + projCost(iter);        
        plotMeNow = plotMe(f) ;
        
        
        t = (sqrt(1+4*tp^2)+1)/2;
        st =  f + (tp-1)/t*(f-fp);
        if norm(f-fp)<1e-8
            break
        end
        fp = f;
        tp = t;
    end
    f_est = f(:);




dataCost(outerIter) = norm(y-(f_est + f_est.*(A*v_est)))^2/2;

fprintf('iter=%d, ||y-yhat||=%e\n',outerIter, dataCost(outerIter))
% plot results
if display && ~mod(outerIter,1)
    figure(100);clf
    subplot(257)
    imagesc(reshape(f_est,imgDim-numRemove,imgDim));colorbar;axis tight equal;
    title(['estimated f'])
    axis off
%     
%     subplot(258)
%     imagesc(reshape(y1,imgDim-numRemove,imgDim));colorbar;axis tight equal;
%     title({'true f';'(Camera shifted for this truth measurment!)'})
% 
%     axis off

    
    subplot(254)
    temp = reshape(P_para*f_est,[imgDim-numRemove,imgDim]);
    imagesc(temp(:,1:end-1))
    colorbar
    colormap('gray')
    axis off
    title('P_{para}*fest')
    
    subplot(253)
    temp = reshape(P_para*y1,[imgDim-numRemove,imgDim]);
    imagesc(temp(:,1:end-1))
    colorbar
    colormap('gray')
    axis off
    title('P_{para}*f')
    
    subplot(2,5,9)
    temp = reshape(P_para*f_est - P_para*y1,[imgDim-numRemove,imgDim]);
    imagesc(temp(:,1:end-1))
    colorbar
    colormap('gray')
    axis off
    title('P_{para}*f_{est} - P_{para}*f')
    
    subplot(2,5,10)
    plot((2:imgDim)/imgDim*90,plotMeAmbient(2:end),'r');
    title('Ambient Light v*')
    
    subplot(255)
    plot((2:imgDim)/imgDim*90,plotMeNow(2:end),'r');
    hold on
    plot((2:imgDim)/imgDim*90,plotMeTruth(2:end),'b')
    xlim([1,90])
    legend('est','truth')
    title({'v*';['normTruth: ' num2str(norm(plotMeTruth(2:end)))]; ['normEst: ' num2str(norm(plotMeNow(2:end)))]})
    
    subplot(251)
    semilogy(dataCost)
    title('Data fidelity cost')
    
    subplot(256)
    semilogy(projCost)
    title('proj cost')
    
    subplot(252)
    plot((1:imgDim)/imgDim*90,sceneWhiteFloor/norm(sceneWhiteFloor),'b.-');
    hold on;
    plot((1:imgDim)/imgDim*90, v_est/norm(v_est),'r.-');
    xlim([1,90])
    legend('reference','estimate','location','best')
    title('estimated v')
    

    drawnow
end
end

  