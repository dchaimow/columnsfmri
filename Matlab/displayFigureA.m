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
set(gca,'FontSize',14);


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
set(gca,'FontSize',14);


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
set(gca,'FontSize',14);
end
