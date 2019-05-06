function results = simulatefMRIOfColumnPatterns(parameters)
% simulatefMRIOfColumnPatterns simulates fMRI of column patterns. 
%   RESULTS = simulatefMRIOfColumnPatterns(PARAMETERS) runs simulations of
%   cortical column patterns and their imaging. It quantifies measures 
%   related to possible analysis objectives and approaches, and estimates
%   optimal voxel widths. Results are displayed as text and returned in 
%   structure RESULTS.
% 
%   PARAMETERS is a structure of parameters that can be set using
%   setParameters(s1,...).
%
%   see also setParameters

% potential limitations of this version:
% - only 2D single slice case
%   alt: isotropic voxels, multiple slices within constant 3D volume
% - no arbitrary spectra
%   alt: specify filter spectrum and/or amplitude
% - detection probability calculated from mvCNR fitting
%   alt: calculate analytically simulate detection probability
% - decoding probability and decoding accuracy from mvCNR fitting
%   alt: simulate decoding
%        NOTE: fitting and simulation depends on some implicit assumptions 
%        (e.g. trial length)
% - output is printout of optimum voxels and associated quantities, in 
%   addition all quantities as a function of voxel width are returned in
%   structure results
%   alt: display of entire curves, display examplary simulated patterns

% set parameters
nTrials        = parameters.nTrials;
Nsim           = parameters.N;
L              = parameters.L;
rho            = parameters.rho;
delta          = parameters.delta;
fwhm           = parameters.fwhm;
b              = parameters.beta;
downFactors    = parameters.downFactors;
sliceThickness = parameters.sliceThickness;
TR             = parameters.TR;
nT             = parameters.nT;
noiseType      = parameters.noiseType;
AFlat          = parameters.AFlat;

% setup simulation grid
sim = setupsim(Nsim,L);

% compute range of voxel widths
wRange = sim.L./(2 .* ceil(sim.N .* downFactors .* 0.5));

% voxel volume as a function of voxel width
voxelVOfW = sliceThickness.*wRange.^2;

% number of voxels as a function of voxel width (assuming single slice)
nVoxelsOfW = AFlat./(wRange.^2);

% noise level as a function of voxel volume
differentialFlag = true;
noiseOfW    = ...
    noiseModel(voxelVOfW,noiseType,TR,nT,differentialFlag);

cr = zeros(nTrials,length(wRange));
cor = zeros(nTrials,length(wRange));

for zTrial = 1:nTrials
    rng(parameters.randomNumberSeed+zTrial);

    % initialize noise patter for simulation of columns
    noise = sim_gwnoise(sim);
    % simulate neuronal columnar pattern
    neuronal = sim_columnPattern(sim,rho,delta,noise); 
    % simulate BOLD response pattern
    bold = sim_bold(sim,fwhm,b,neuronal);
    
    % for each tested voxel width
    for zW = 1:length(wRange)
        w = wRange(zW);
        
        % simulate MRI voxel sampling
        downFactor = sim.dx/w;
        [voxel,~] = sim_mri(sim,downFactor,bold);
        
        % calculate contrast range
        cr(zTrial,zW) = std(voxel(:));
        
        % add noise
        voxelPlusNoise = voxel + noiseOfW(zW)*randn(size(voxel));
   
        % calculate correlation to original (neuronal) pattern        
        cor(zTrial,zW) =  patternCorrelation(sim,neuronal,voxelPlusNoise);
    end        
end

% average over trials
cr = mean(cr);
cor = mean(cor);

% calculate CNR and mv-CNR
CNR = cr./noiseOfW;
a = 0.26;
mvCNR = CNR.*nVoxelsOfW.^a;

% calculate univariate detection probability
pUniOfW = detectionProbability(CNR,1);

% calculate multivariate detection probability
fPMultiDetect = cdfModel(1.9038,3.5990,0.05,@wblcdf);
pMultiOfW = fPMultiDetect(mvCNR);

% calculate decoding probability
fPDecode2 = cdfModel(2.0675,3.4101,0.05,@wblcdf);
fPDecode4 = cdfModel(2.8451,3.4770,0.05,@wblcdf);
fPDecode8 = cdfModel(3.8086,3.2733,0.05,@wblcdf);

pDecodeOfW2 = fPDecode2(mvCNR);
pDecodeOfW4 = fPDecode4(mvCNR);
pDecodeOfW8 = fPDecode8(mvCNR);

% calculate decoding accuracy
fMeanClassPerf2 = cdfModel(4.7018,1.9954,1/2,@wblcdf);
fMeanClassPerf4 = cdfModel(8.0900,2.1585,1/4,@wblcdf);
fMeanClassPerf8 = cdfModel(15.8022,1.9954,1/8,@wblcdf);

meanClassPerfOfW2 = fMeanClassPerf2(mvCNR);
meanClassPerfOfW4 = fMeanClassPerf4(mvCNR);
meanClassPerfOfW8 = fMeanClassPerf8(mvCNR);

% find optima
[maxPUni,idxMaxPUni] = max(pUniOfW);
optWPUni = wRange(idxMaxPUni);

[maxMVCNR,idxMaxMVCNR] = max(mvCNR);
optWMVCNR = wRange(idxMaxMVCNR);

maxPMulti = max(pMultiOfW);
assert(maxPMulti==pMultiOfW(idxMaxMVCNR));

maxPDecode2 = max(pDecodeOfW2);
maxPDecode4 = max(pDecodeOfW4);
maxPDecode8 = max(pDecodeOfW8);
assert(maxPDecode2==pDecodeOfW2(idxMaxMVCNR));
assert(maxPDecode4==pDecodeOfW4(idxMaxMVCNR));
assert(maxPDecode8==pDecodeOfW8(idxMaxMVCNR));

maxMeanClassPerf2 = max(meanClassPerfOfW2);
maxMeanClassPerf4 = max(meanClassPerfOfW4);
maxMeanClassPerf8 = max(meanClassPerfOfW8);
assert(maxMeanClassPerf2==meanClassPerfOfW2(idxMaxMVCNR));
assert(maxMeanClassPerf4==meanClassPerfOfW4(idxMaxMVCNR));
assert(maxMeanClassPerf8==meanClassPerfOfW8(idxMaxMVCNR));

[maxCor,idxMaxCor] = max(cor);
optWCor = wRange(idxMaxCor);

% assign results
results.w                        = wRange;
results.cr                       = cr;
results.noiseLevel               = noiseOfW;
results.SNR                      = 1./noiseOfW;
results.CNR                      = CNR;

results.pDetectUnivariate        = pUniOfW;
results.pDetectUnivariate_max    = maxPUni;
results.pDetectUnivariate_optW   = optWPUni;

results.mvCNR                    = mvCNR;
results.mvCNR_max                = maxMVCNR;
results.mvCNR_opt                = optWMVCNR;

results.pDetectMultivariate      = pMultiOfW;
results.pDetectMultivariate_max  = maxPMulti;
results.pDetectMultivariate_optW = optWMVCNR;

results.pDecode2                  = pDecodeOfW2;
results.pDecode2_max              = maxPDecode2;
results.pDecode2_optW             = optWMVCNR;

results.accuracyDecode2           = meanClassPerfOfW2;
results.accuracyDecode2_max       = maxMeanClassPerf2;
results.accuracyDecode2_optW      = optWMVCNR;

results.pDecode4                  = pDecodeOfW4;
results.pDecode4_max              = maxPDecode4;
results.pDecode4_optW             = optWMVCNR;

results.accuracyDecode4           = meanClassPerfOfW4;
results.accuracyDecode4_max       = maxMeanClassPerf4;
results.accuracyDecode4_optW      = optWMVCNR;

results.pDecode8                  = pDecodeOfW8;
results.pDecode8_max              = maxPDecode8;
results.pDecode8_optW             = optWMVCNR;

results.accuracyDecode8           = meanClassPerfOfW8;
results.accuracyDecode8_max       = maxMeanClassPerf8;
results.accuracyDecode8_optW      = optWMVCNR;

results.patternCorrelation       = cor;
results.patternCorrelation_max   = maxCor;
results.patternCorrelation_optW  = optWCor;

% print and display results
printResults(results);
displayFigureA(results);
displayFigureB(results);
displayFigureC(results);
end

function printResults(results)
% print out results
fprintf('optimized quantity                 | optimal value | optimal voxel width\n');
fprintf('-----------------------------------+---------------+--------------------\n');
fprintf('univariate detection probability   | %.2f          | %.2f mm            \n',...
    results.pDetectUnivariate_max,results.pDetectUnivariate_optW);
fprintf('multivariate detection probability | %.2f          | %.2f mm            \n',...
    results.pDetectMultivariate_max,results.pDetectMultivariate_optW);
fprintf('decoding probability - 2 classes   | %.2f          | %.2f mm            \n',...
    results.pDecode2_max,results.pDecode2_optW);
fprintf('decoding accuracy    - 2 classes   | %.2f          | %.2f mm            \n',...
    results.accuracyDecode2_max,results.accuracyDecode2_optW);
fprintf('decoding probability - 4 classes   | %.2f          | %.2f mm            \n',...
    results.pDecode4_max,results.pDecode4_optW);
fprintf('decoding accuracy    - 4 classes   | %.2f          | %.2f mm            \n',...
    results.accuracyDecode4_max,results.accuracyDecode4_optW);
fprintf('decoding probability - 8 classes   | %.2f          | %.2f mm            \n',...
    results.pDecode8_max,results.pDecode8_optW);
fprintf('decoding accuracy    - 8 classes   | %.2f          | %.2f mm            \n',...
    results.accuracyDecode8_max,results.accuracyDecode8_optW);
fprintf('pattern correlation                | %.2f          | %.2f mm            \n',...
    results.patternCorrelation_max,results.patternCorrelation_optW);
end

function displayFigureA(results)
% figure A - optimization of all 3 objectives
f = createFigure();
c = get(f,'DefaultAxesColorOrder');

subplot(1,3,1);
plot(results.w,results.pDetectUnivariate,'Color',c(2,:));
line([results.pDetectUnivariate_optW results.pDetectUnivariate_optW],...
    [0 1],'Color',c(2,:),'LineStyle','--')
axis([min(results.w) max(results.w) 0 1]);
axis square;
box off;
xlabel('voxel width [mm]');
ylabel('probability');
title('Univariate detection');

subplot(1,3,2);
plot(results.w,results.accuracyDecode2,'Color',c(5,:));
line([results.accuracyDecode2_optW results.accuracyDecode2_optW],[0 1],...
    'Color',c(5,:),'LineStyle','--')
axis([min(results.w) max(results.w) 0.5 1]);
axis square;
box off;
xlabel('voxel width [mm]');
ylabel('accuracy');
title('Multivariate decoding');

subplot(1,3,3);
plot(results.w,results.patternCorrelation,'Color',c(1,:));
line([results.patternCorrelation_optW results.patternCorrelation_optW],...
    [0 1],'Color',c(1,:),'LineStyle','--')
axis([min(results.w) max(results.w) 0 1]);
axis square;
box off;
xlabel('voxel width [mm]');
ylabel('correlation');
title('Reconstruction');
end

function displayFigureB(results)
% figure B - optimization of single voxel CNR
f = createFigure();
c = get(f,'DefaultAxesColorOrder');

subplot(1,2,1);
ax = visualize1d2func(results.w,[100*results.cr; results.SNR]',...
    'voxel width [mm]',...
    {'signal change [%]','SNR'},'',{'k',c(3,:)});
axis(ax(1),[0 4 0 max(100*results.cr)*1.5],'square');
axis(ax(2),[0 4 0 max(results.SNR)*1.5],'square');
set(ax(1),'YTick',[0 0.5 1 1.5 2]);
set(ax(2),'YTick',[0 100 200 300 400 500]);
box off;
title(...
    {'Contrast range and distribution of '...
    'differential responses vs. SNR'});

subplot(1,2,2);
nFractions = 5;
plotNormPFractions(results.w,results.CNR',nFractions,hsv2rgb([0.1 1 0.5]));
hold on;
xlabel('voxel width [mm]');
ylabel('response ratio');

axis([0 4 0 1.5*max(results.CNR)],'square');
box off;
title({'CNR and distribution of ',...
    'differential responses relative to noise'});
[~,maxCNRIdx] = max(results.CNR);
line(results.w(maxCNRIdx)*[1 1],[0 1.5*max(results.CNR)],...
    'Color',hsv2rgb([0.1 1 0.5]),'LineStyle',':');
end

function displayFigureC(results)
% figure C - optimization of multivariate analysis
f = createFigure();
c = get(f,'DefaultAxesColorOrder');

42;

end

function displayFigureD(results)
% figure D - optimization of reconstruction
setFigureOptions;
c = get(0,'DefaultAxesColorOrder');

end

function f = createFigure(figureSize)
c = [...
    0.6   1      0.7    ;...
    1      1  0.9   ;...
    0.1  1 1      ;...
    0.8 1 0.6    ;...
    0.3333  1 0.5;...
    0.58   0.67   0.9];

screenSize = get(0,'ScreenSize');

f = figure;

% set line properties:
set(f,'DefaultLineLineWidth',1.5);
set(f,'DefaultLineMarkerSize',4);

% set axes Properties
set(f,'DefaultAxesLineWidth',1);
set(f,'DefaultAxesFontSize',12);
set(f,'DefaultAxesFontName','Arial');
set(f,'DefaultAxesColor','None');
set(f,'DefaultAxesColorOrder',hsv2rgb(c));

% set figure properties
set(f,'Color','White');
figurePos = get(f,'Position');
set(f,'Position',...
    [figurePos(1) screenSize(4)-figurePos(2),...
    figurePos(3:4)]);
end

function f = cdfModel(a,b,y0,cdf)
f = @(x) y0+(cdf(x,a,b)-cdf(0,a,b))*(1-y0);
end

function AX=visualize1d2func(xRange,data,xLabel,yLabel,titleText,lineStyle)
nFractions = 5;
color = [0 0 0];
plotFunction1 = @(x,y) plotNormPFractions(x,y,nFractions,color);
plotFunction2 = 'plot';

[AX,H1,H2] = plotyy(xRange, data(:,1),xRange,data(:,2),plotFunction1,plotFunction2);
hold on;
set(get(AX(1),'Ylabel'),'String',yLabel(1));
%set(get(AX(1),'Ylabel'));
set(get(AX(2),'Ylabel'),'String',yLabel(2));
%set(get(AX(2),'Ylabel'));
xlabel(xLabel);
if exist('lineStyle','var')
%    set(H1,'Color',lineStyle{1});
%    set(AX(1),'YColor',lineStyle{1});
    set(get(AX(1),'Ylabel'),'Color',lineStyle{1});
    set(H2,'Color',lineStyle{2});
    set(AX(2),'YColor',lineStyle{2})
    set(get(AX(2),'Ylabel'),'Color',lineStyle{2});
end
title(titleText);
uistack(AX(2));
end

function visualize1d(xRange,data,xLabel,yLabel,titleText,lineStyle)
if ~exist('lineStyle','var')
    lineStyle = 'k';
end

if isnumeric(lineStyle)
    plot(xRange, data,'Color',lineStyle);
else
    plot(xRange, data,lineStyle);
end
xlabel(xLabel);
ylabel(yLabel);
title(titleText);
end

function putGridLines(xRes,yRes,leaveOutLast)
ax = axis(gca);
if ~isempty(xRes)
    startX = floor(ax(1)/xRes)*xRes;
    endX = ceil(ax(2)/xRes)*xRes;
    if exist('leaveOutLast','var') && leaveOutLast
        endX = endX - xRes;
    end
    for x=startX:xRes:endX
        l = line([x x],[ax(3) ax(4)],'LineStyle',':','Color','k','LineWidth',1);
        uistack(l,'bottom');
    end
end
if ~isempty(yRes)
    startY = ceil(ax(3)/yRes)*yRes;
    endY = ceil(ax(4)/yRes)*yRes;
    for y=startY:yRes:endY
        l = line([ax(1) ax(2)],[y y],'LineStyle',':','Color','k','LineWidth',1);
        uistack(l,'bottom');
    end
end
end

function h = plotNormPFractions(x,y,nFractions,c)
yFractions = norminv(0.5+ (0:1/nFractions:1)/2,0,1);
cGray = (0.6*c + 0.4*[1 1 1]);
cWhite = [1 1 1];
for z=1:nFractions
    a = (z-1)/(nFractions);
    cArea = cWhite * a + (1 - a) * cGray;
    plot_area(x,...
        y'*yFractions(z+1),y'*yFractions(z),cArea);
    hold on;
end
h = plot(x,y,'Color',c);
hold off;
end

function h = plot_area(x,lower,upper,color)
lower(isinf(lower)) = 10*sign(lower(isinf(lower))); % check if 10 factor is large enough
upper(isinf(upper)) = 10*sign(lower(isinf(upper))); % check if 10 factor is large enough

h = set(fill([x,x(end:-1:1)],...
    [upper,lower(end:-1:1)],color),'EdgeColor',color);
end