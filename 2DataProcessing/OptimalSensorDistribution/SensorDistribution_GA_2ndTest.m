clear 
FS = 15;
colors = lines(7);

NumSen=8; % Set the number of selected sensors
NumSen2 = NumSen-1;

%% =========================== Load results ===============================

% addpath('add the location of the following files')  
addpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/2DataProcessing/FIM-FIMADPR-Assessment')
addpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/2DataProcessing/SSC-SSCADPR-Assessment')

load FIMGASensorDistribution.mat bestIndicesFIMGA DelNumGA bestIndicesFIMFSGA2 DelNumGAFS bestIndicesFIMFSRGA2 DelNumGAFSR RepNumGAFSR
DelNumFIMGA=DelNumGA; DelNumFIMGAFS=DelNumGAFS; 
DelNumFIMGAFSR=DelNumGAFSR; RepNumFIMGAFSR=RepNumGAFSR;
clear DelNumGA DelNumGAFS DelNumGAFSR RepNumGAFSR

load FIMADPRGASensorDistribution.mat bestIndicesFIMADPRGA DelNumADPRGA bestIndicesFIMADPR_FSGA2 DelNumADPRGAFS bestIndicesFIMADPR_FSRGA2 DelNumADPRGAFSR RepNumADPRGAFSR
DelNumFIMADPRGA=DelNumADPRGA; DelNumFIMADPRGAFS=DelNumADPRGAFS;
DelNumFIMADPRGAFSR=DelNumADPRGAFSR;RepNumFIMADPRGAFSR=RepNumADPRGAFSR;
clear DelNumADPRGA DelNumADPRGAFS DelNumADPRGAFSR RepNumADPRGAFSR

load SSCGASensorDistribution.mat bestIndicesSSCGA DelNumGA bestIndicesSSCFSGA2 DelNumGAFS bestIndicesSSCFSRGA2 DelNumGAFSR RepNumGAFSR
DelNumSSCGA=DelNumGA; DelNumSSCGAFS=DelNumGAFS; 
DelNumSSCGAFSR=DelNumGAFSR; RepNumSSCGAFSR=RepNumGAFSR;
clear DelNumGA DelNumGAFS DelNumGAFSR RepNumGAFSR

load SSCADPRGASensorDistribution.mat bestIndicesSSCADPRGA DelNumGA bestIndicesSSCADPRFSGA2 DelNumGAFS bestIndicesSSCADPRFSRGA2 DelNumGAFSR RepNumGAFSR
DelNumSSCADPRGA=DelNumGA; DelNumSSCADPRGAFS=DelNumGAFS; 
DelNumSSCADPRGAFSR=DelNumGAFSR; RepNumSSCADPRGAFSR=RepNumGAFSR;
clear DelNumGA DelNumGAFS DelNumGAFSR RepNumGAFSR

% rmpath('add the location of these files')
rmpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/2DataProcessing/FIM-FIMADPR-Assessment')
rmpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/2DataProcessing/SSC-SSCADPR-Assessment')

%% ======== Extract the coordinates of sensors ============================
% addpath('add the location of the folder ModeShapeData')  
addpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/ModeShapeData')
load ModeShapeF_T15UD2.mat

% DOF coordinates in space
x= modal_data.Geometry.XYZ(1:36,2);
y =modal_data.Geometry.XYZ(1:36,1);

% rmpath('add the location of the folder ModeShapeData')
rmpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/ModeShapeData')

%% ============================= FIM ======================================

sf = 180;
gap = [0 0]; % [vertical,horizontal]
marg_h = [0.06 0.05]; % [lower uppper]
marg_w = [0.02 0.02]; % [left right]
subplot = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
figure('rend','painters','pos',[100 100 sf*3.4 sf])

% the outline of the wing 
figure(1)
subplot(1, 1, 1)
hold on
XC(1,1) = x(1,1)-0.15; XC(2,1) = x(12,1)+0.15;
YC(1,1) = y(1,1)-0.25; YC(2,1) = y(12,1)-0.25;
plot(XC([1 2],1),YC([1 2],1),'k-','LineWidth',1);

XC(3,1) = x(25,1)-0.15; XC(4,1) = x(31,1); XC(5,1) = x(32,1); XC(6,1) = x(36,1)+0.15;
YC(3,1) = y(25,1)+0.25; YC(4,1) = y(31,1)+0.25; YC(5,1) = y(32,1)+0.25; YC(6,1) = y(36,1)+0.25;

plot(XC(3:6,1),YC(3:6,1),'k-','LineWidth',1)
plot(XC([1 3],1),YC([1 3],1),'k-','LineWidth',1)
plot(XC([2 6],1),YC([2 6],1),'k-','LineWidth',1)

% Sensor location coordinates and labels

plot(x,y,'ks','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');

% =========================================================================
plot(x,y,'rd','MarkerIndices',[bestIndicesFIMGA{NumSen,1}],'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7);  
plot(x,y,'rd','MarkerIndices',bestIndicesFIMGA{NumSen,1}(1,DelNumFIMGA(1,NumSen)),'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerSize',7); 

set(gca,'XTick',[], 'YTick', [])
axis equal
axis([-0.5 8 -0.5 1.5])
hold off

%% ============================= FIM-FS ===================================

sf = 180;
gap = [0 0]; % [vertical,horizontal]
marg_h = [0.06 0.05]; % [lower uppper]
marg_w = [0.02 0.02]; % [left right]
subplot = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
figure('rend','painters','pos',[100 100 sf*3.4 sf])

% the outline of the wing 
figure(2)
subplot(1, 1, 1)
hold on
XC(1,1) = x(1,1)-0.15; XC(2,1) = x(12,1)+0.15;
YC(1,1) = y(1,1)-0.25; YC(2,1) = y(12,1)-0.25;
plot(XC([1 2],1),YC([1 2],1),'k-','LineWidth',1);

XC(3,1) = x(25,1)-0.15; XC(4,1) = x(31,1); XC(5,1) = x(32,1); XC(6,1) = x(36,1)+0.15;
YC(3,1) = y(25,1)+0.25; YC(4,1) = y(31,1)+0.25; YC(5,1) = y(32,1)+0.25; YC(6,1) = y(36,1)+0.25;

plot(XC(3:6,1),YC(3:6,1),'k-','LineWidth',1)
plot(XC([1 3],1),YC([1 3],1),'k-','LineWidth',1)
plot(XC([2 6],1),YC([2 6],1),'k-','LineWidth',1)

% Sensor location coordinates and labels

plot(x,y,'ks','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');

% =========================================================================
plot(x,y,'rd','MarkerIndices',[bestIndicesFIMFSGA2{NumSen,1}],'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7);  
plot(x,y,'rd','MarkerIndices',bestIndicesFIMFSGA2{NumSen,1}(1,DelNumFIMGAFS(1,NumSen)),'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerSize',7); 

set(gca,'XTick',[], 'YTick', [])
axis equal
axis([-0.5 8 -0.5 1.5])
hold off

%% ============================= FIM-FSR ==================================

sf = 180;
gap = [0 0]; % [vertical,horizontal]
marg_h = [0.06 0.05]; % [lower uppper]
marg_w = [0.02 0.02]; % [left right]
subplot = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
figure('rend','painters','pos',[100 100 sf*3.4 sf])

% the outline of the wing 
figure(3)
subplot(1, 1, 1)
hold on
XC(1,1) = x(1,1)-0.15; XC(2,1) = x(12,1)+0.15;
YC(1,1) = y(1,1)-0.25; YC(2,1) = y(12,1)-0.25;
plot(XC([1 2],1),YC([1 2],1),'k-','LineWidth',1);

XC(3,1) = x(25,1)-0.15; XC(4,1) = x(31,1); XC(5,1) = x(32,1); XC(6,1) = x(36,1)+0.15;
YC(3,1) = y(25,1)+0.25; YC(4,1) = y(31,1)+0.25; YC(5,1) = y(32,1)+0.25; YC(6,1) = y(36,1)+0.25;

plot(XC(3:6,1),YC(3:6,1),'k-','LineWidth',1)
plot(XC([1 3],1),YC([1 3],1),'k-','LineWidth',1)
plot(XC([2 6],1),YC([2 6],1),'k-','LineWidth',1)

% Sensor location coordinates and labels

plot(x,y,'ks','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');

% =========================================================================
plot(x,y,'rd','MarkerIndices',[bestIndicesFIMFSRGA2{NumSen2,1}],'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7);  
plot(x,y,'rd','MarkerIndices',bestIndicesFIMFSRGA2{NumSen2,1}(1,DelNumFIMGAFSR(1,NumSen2)),'MarkerEdgeColor',colors(5,:),'MarkerFaceColor',colors(5,:),'MarkerSize',7); 
plot(x+0.1,y,'rd','MarkerIndices',bestIndicesFIMFSRGA2{NumSen2,1}(1,RepNumFIMGAFSR(1,NumSen2)),'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7); 

set(gca,'XTick',[], 'YTick', [])
axis equal
axis([-0.5 8 -0.5 1.5])
hold off

%% ============================= FIM-ADPR =================================

sf = 180;
gap = [0 0]; % [vertical,horizontal]
marg_h = [0.06 0.05]; % [lower uppper]
marg_w = [0.02 0.02]; % [left right]
subplot = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
figure('rend','painters','pos',[100 100 sf*3.4 sf])

% the outline of the wing 
figure(4)
subplot(1, 1, 1)
hold on
XC(1,1) = x(1,1)-0.15; XC(2,1) = x(12,1)+0.15;
YC(1,1) = y(1,1)-0.25; YC(2,1) = y(12,1)-0.25;
plot(XC([1 2],1),YC([1 2],1),'k-','LineWidth',1);

XC(3,1) = x(25,1)-0.15; XC(4,1) = x(31,1); XC(5,1) = x(32,1); XC(6,1) = x(36,1)+0.15;
YC(3,1) = y(25,1)+0.25; YC(4,1) = y(31,1)+0.25; YC(5,1) = y(32,1)+0.25; YC(6,1) = y(36,1)+0.25;

plot(XC(3:6,1),YC(3:6,1),'k-','LineWidth',1)
plot(XC([1 3],1),YC([1 3],1),'k-','LineWidth',1)
plot(XC([2 6],1),YC([2 6],1),'k-','LineWidth',1)

% Sensor location coordinates and labels

plot(x,y,'ks','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');

% =========================================================================
plot(x,y,'rd','MarkerIndices',[bestIndicesFIMADPRGA{NumSen,1}],'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7);  
plot(x,y,'rd','MarkerIndices',bestIndicesFIMADPRGA{NumSen,1}(1,DelNumFIMADPRGA(1,NumSen)),'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerSize',7); 

set(gca,'XTick',[], 'YTick', [])
axis equal
axis([-0.5 8 -0.5 1.5])
hold off

%% ============================= FIM-ADPR-FS ==============================

sf = 180;
gap = [0 0]; % [vertical,horizontal]
marg_h = [0.06 0.05]; % [lower uppper]
marg_w = [0.02 0.02]; % [left right]
subplot = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
figure('rend','painters','pos',[100 100 sf*3.4 sf])

% the outline of the wing 
figure(5)
subplot(1, 1, 1)
hold on
XC(1,1) = x(1,1)-0.15; XC(2,1) = x(12,1)+0.15;
YC(1,1) = y(1,1)-0.25; YC(2,1) = y(12,1)-0.25;
plot(XC([1 2],1),YC([1 2],1),'k-','LineWidth',1);

XC(3,1) = x(25,1)-0.15; XC(4,1) = x(31,1); XC(5,1) = x(32,1); XC(6,1) = x(36,1)+0.15;
YC(3,1) = y(25,1)+0.25; YC(4,1) = y(31,1)+0.25; YC(5,1) = y(32,1)+0.25; YC(6,1) = y(36,1)+0.25;

plot(XC(3:6,1),YC(3:6,1),'k-','LineWidth',1)
plot(XC([1 3],1),YC([1 3],1),'k-','LineWidth',1)
plot(XC([2 6],1),YC([2 6],1),'k-','LineWidth',1)

% Sensor location coordinates and labels

plot(x,y,'ks','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');

% =========================================================================
plot(x,y,'rd','MarkerIndices',[bestIndicesFIMADPR_FSGA2{NumSen,1}],'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7);  
plot(x,y,'rd','MarkerIndices',bestIndicesFIMADPR_FSGA2{NumSen,1}(1,DelNumFIMADPRGAFS(1,NumSen)),'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerSize',7); 

set(gca,'XTick',[], 'YTick', [])
axis equal
axis([-0.5 8 -0.5 1.5])
hold off

%% ============================= FIM-ADPR-FSR =============================

sf = 180;
gap = [0 0]; % [vertical,horizontal]
marg_h = [0.06 0.05]; % [lower uppper]
marg_w = [0.02 0.02]; % [left right]
subplot = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
figure('rend','painters','pos',[100 100 sf*3.4 sf])

% the outline of the wing 
figure(6)
subplot(1, 1, 1)
hold on
XC(1,1) = x(1,1)-0.15; XC(2,1) = x(12,1)+0.15;
YC(1,1) = y(1,1)-0.25; YC(2,1) = y(12,1)-0.25;
plot(XC([1 2],1),YC([1 2],1),'k-','LineWidth',1);

XC(3,1) = x(25,1)-0.15; XC(4,1) = x(31,1); XC(5,1) = x(32,1); XC(6,1) = x(36,1)+0.15;
YC(3,1) = y(25,1)+0.25; YC(4,1) = y(31,1)+0.25; YC(5,1) = y(32,1)+0.25; YC(6,1) = y(36,1)+0.25;

plot(XC(3:6,1),YC(3:6,1),'k-','LineWidth',1)
plot(XC([1 3],1),YC([1 3],1),'k-','LineWidth',1)
plot(XC([2 6],1),YC([2 6],1),'k-','LineWidth',1)

% Sensor location coordinates and labels

plot(x,y,'ks','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');

% =========================================================================
plot(x,y,'rd','MarkerIndices',[bestIndicesFIMADPR_FSRGA2{NumSen2,1}],'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7);  
plot(x,y,'rd','MarkerIndices',bestIndicesFIMADPR_FSRGA2{NumSen2,1}(1,DelNumFIMADPRGAFSR(1,NumSen2)),'MarkerEdgeColor',colors(5,:),'MarkerFaceColor',colors(5,:),'MarkerSize',7); 
plot(x+0.1,y,'rd','MarkerIndices',bestIndicesFIMADPR_FSRGA2{NumSen2,1}(1,RepNumFIMADPRGAFSR(1,NumSen2)),'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7);  

set(gca,'XTick',[], 'YTick', [])
axis equal
axis([-0.5 8 -0.5 1.5])
hold off

%% ============================= SSC ======================================

sf = 180;
gap = [0 0]; % [vertical,horizontal]
marg_h = [0.06 0.05]; % [lower uppper]
marg_w = [0.02 0.02]; % [left right]
subplot = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
figure('rend','painters','pos',[100 100 sf*3.4 sf])

% the outline of the wing 
figure(7)
subplot(1, 1, 1)
hold on
XC(1,1) = x(1,1)-0.15; XC(2,1) = x(12,1)+0.15;
YC(1,1) = y(1,1)-0.25; YC(2,1) = y(12,1)-0.25;
plot(XC([1 2],1),YC([1 2],1),'k-','LineWidth',1);

XC(3,1) = x(25,1)-0.15; XC(4,1) = x(31,1); XC(5,1) = x(32,1); XC(6,1) = x(36,1)+0.15;
YC(3,1) = y(25,1)+0.25; YC(4,1) = y(31,1)+0.25; YC(5,1) = y(32,1)+0.25; YC(6,1) = y(36,1)+0.25;

plot(XC(3:6,1),YC(3:6,1),'k-','LineWidth',1)
plot(XC([1 3],1),YC([1 3],1),'k-','LineWidth',1)
plot(XC([2 6],1),YC([2 6],1),'k-','LineWidth',1)

% Sensor location coordinates and labels
plot(x,y,'ks','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');

% Damage location coordinates and labels
M1x=(x(4,1)+x(5,1))/2; M1y=(y(4,1)+y(5,1))/2;
M2x=(x(6,1)+x(7,1))/2; M2y=(y(6,1)+y(7,1))/2;
M3x=(x(11,1)+x(12,1))/2; M3y=(y(11,1)+y(12,1))/2;
plot(M1x+0.05,M1y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
plot(M2x+0.05,M2y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
plot(M3x+0.05,M3y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);

% =========================================================================
plot(x,y,'rd','MarkerIndices',[bestIndicesSSCGA{NumSen,1}],'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7);  
plot(x,y,'rd','MarkerIndices',bestIndicesSSCGA{NumSen,1}(1,DelNumSSCGA(1,NumSen)),'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerSize',7); 

set(gca,'XTick',[], 'YTick', [])
axis equal
axis([-0.5 8 -0.5 1.5])
hold off

%% ============================= SSC-FS ===================================

sf = 180;
gap = [0 0]; % [vertical,horizontal]
marg_h = [0.06 0.05]; % [lower uppper]
marg_w = [0.02 0.02]; % [left right]
subplot = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
figure('rend','painters','pos',[100 100 sf*3.4 sf])

% the outline of the wing 
figure(8)
subplot(1, 1, 1)
hold on
XC(1,1) = x(1,1)-0.15; XC(2,1) = x(12,1)+0.15;
YC(1,1) = y(1,1)-0.25; YC(2,1) = y(12,1)-0.25;
plot(XC([1 2],1),YC([1 2],1),'k-','LineWidth',1);

XC(3,1) = x(25,1)-0.15; XC(4,1) = x(31,1); XC(5,1) = x(32,1); XC(6,1) = x(36,1)+0.15;
YC(3,1) = y(25,1)+0.25; YC(4,1) = y(31,1)+0.25; YC(5,1) = y(32,1)+0.25; YC(6,1) = y(36,1)+0.25;

plot(XC(3:6,1),YC(3:6,1),'k-','LineWidth',1)
plot(XC([1 3],1),YC([1 3],1),'k-','LineWidth',1)
plot(XC([2 6],1),YC([2 6],1),'k-','LineWidth',1)

% Sensor location coordinates and labels
plot(x,y,'ks','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');

% Damage location coordinates and labels
M1x=(x(4,1)+x(5,1))/2; M1y=(y(4,1)+y(5,1))/2;
M2x=(x(6,1)+x(7,1))/2; M2y=(y(6,1)+y(7,1))/2;
M3x=(x(11,1)+x(12,1))/2; M3y=(y(11,1)+y(12,1))/2;
plot(M1x+0.05,M1y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
plot(M2x+0.05,M2y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
plot(M3x+0.05,M3y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);

% =========================================================================
plot(x,y,'rd','MarkerIndices',[bestIndicesSSCFSGA2{NumSen,1}],'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7);  
plot(x,y,'rd','MarkerIndices',bestIndicesSSCFSGA2{NumSen,1}(1,DelNumSSCGAFS(1,NumSen)),'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerSize',7); 

set(gca,'XTick',[], 'YTick', [])
axis equal
axis([-0.5 8 -0.5 1.5])
hold off

%% ============================= SSC-FSR ==================================

sf = 180;
gap = [0 0]; % [vertical,horizontal]
marg_h = [0.06 0.05]; % [lower uppper]
marg_w = [0.02 0.02]; % [left right]
subplot = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
figure('rend','painters','pos',[100 100 sf*3.4 sf])

% the outline of the wing 
figure(9)
subplot(1, 1, 1)
hold on
XC(1,1) = x(1,1)-0.15; XC(2,1) = x(12,1)+0.15;
YC(1,1) = y(1,1)-0.25; YC(2,1) = y(12,1)-0.25;
plot(XC([1 2],1),YC([1 2],1),'k-','LineWidth',1);

XC(3,1) = x(25,1)-0.15; XC(4,1) = x(31,1); XC(5,1) = x(32,1); XC(6,1) = x(36,1)+0.15;
YC(3,1) = y(25,1)+0.25; YC(4,1) = y(31,1)+0.25; YC(5,1) = y(32,1)+0.25; YC(6,1) = y(36,1)+0.25;

plot(XC(3:6,1),YC(3:6,1),'k-','LineWidth',1)
plot(XC([1 3],1),YC([1 3],1),'k-','LineWidth',1)
plot(XC([2 6],1),YC([2 6],1),'k-','LineWidth',1)

% Sensor location coordinates and labels
plot(x,y,'ks','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');

% Damage location coordinates and labels
M1x=(x(4,1)+x(5,1))/2; M1y=(y(4,1)+y(5,1))/2;
M2x=(x(6,1)+x(7,1))/2; M2y=(y(6,1)+y(7,1))/2;
M3x=(x(11,1)+x(12,1))/2; M3y=(y(11,1)+y(12,1))/2;
plot(M1x+0.05,M1y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
plot(M2x+0.05,M2y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
plot(M3x+0.05,M3y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);

% =========================================================================
plot(x,y,'rd','MarkerIndices',[bestIndicesSSCFSRGA2{NumSen2,1}],'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7);  
plot(x,y,'rd','MarkerIndices',bestIndicesSSCFSRGA2{NumSen2,1}(1,DelNumSSCGAFSR(1,NumSen2)),'MarkerEdgeColor',colors(5,:),'MarkerFaceColor',colors(5,:),'MarkerSize',7); 
plot(x+0.1,y,'rd','MarkerIndices',bestIndicesSSCFSRGA2{NumSen2,1}(1,RepNumSSCGAFSR(1,NumSen2)),'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7); 

set(gca,'XTick',[], 'YTick', [])
axis equal
axis([-0.5 8 -0.5 1.5])
hold off

%% ============================= SSC-ADPR =================================

sf = 180;
gap = [0 0]; % [vertical,horizontal]
marg_h = [0.06 0.05]; % [lower uppper]
marg_w = [0.02 0.02]; % [left right]
subplot = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
figure('rend','painters','pos',[100 100 sf*3.4 sf])

% the outline of the wing 
figure(10)
subplot(1, 1, 1)
hold on
XC(1,1) = x(1,1)-0.15; XC(2,1) = x(12,1)+0.15;
YC(1,1) = y(1,1)-0.25; YC(2,1) = y(12,1)-0.25;
plot(XC([1 2],1),YC([1 2],1),'k-','LineWidth',1);

XC(3,1) = x(25,1)-0.15; XC(4,1) = x(31,1); XC(5,1) = x(32,1); XC(6,1) = x(36,1)+0.15;
YC(3,1) = y(25,1)+0.25; YC(4,1) = y(31,1)+0.25; YC(5,1) = y(32,1)+0.25; YC(6,1) = y(36,1)+0.25;

plot(XC(3:6,1),YC(3:6,1),'k-','LineWidth',1)
plot(XC([1 3],1),YC([1 3],1),'k-','LineWidth',1)
plot(XC([2 6],1),YC([2 6],1),'k-','LineWidth',1)

% Sensor location coordinates and labels
plot(x,y,'ks','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');

% Damage location coordinates and labels
M1x=(x(4,1)+x(5,1))/2; M1y=(y(4,1)+y(5,1))/2;
M2x=(x(6,1)+x(7,1))/2; M2y=(y(6,1)+y(7,1))/2;
M3x=(x(11,1)+x(12,1))/2; M3y=(y(11,1)+y(12,1))/2;
plot(M1x+0.05,M1y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
plot(M2x+0.05,M2y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
plot(M3x+0.05,M3y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);

% =========================================================================
plot(x,y,'rd','MarkerIndices',[bestIndicesSSCADPRGA{NumSen,1}],'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7);  
plot(x,y,'rd','MarkerIndices',bestIndicesSSCADPRGA{NumSen,1}(1,DelNumSSCADPRGA(1,NumSen)),'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerSize',7); 

set(gca,'XTick',[], 'YTick', [])
axis equal
axis([-0.5 8 -0.5 1.5])
hold off

%% ============================= SSC-ADPR-FS ==============================

sf = 180;
gap = [0 0]; % [vertical,horizontal]
marg_h = [0.06 0.05]; % [lower uppper]
marg_w = [0.02 0.02]; % [left right]
subplot = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
figure('rend','painters','pos',[100 100 sf*3.4 sf])

% the outline of the wing 
figure(11)
subplot(1, 1, 1)
hold on
XC(1,1) = x(1,1)-0.15; XC(2,1) = x(12,1)+0.15;
YC(1,1) = y(1,1)-0.25; YC(2,1) = y(12,1)-0.25;
plot(XC([1 2],1),YC([1 2],1),'k-','LineWidth',1);

XC(3,1) = x(25,1)-0.15; XC(4,1) = x(31,1); XC(5,1) = x(32,1); XC(6,1) = x(36,1)+0.15;
YC(3,1) = y(25,1)+0.25; YC(4,1) = y(31,1)+0.25; YC(5,1) = y(32,1)+0.25; YC(6,1) = y(36,1)+0.25;

plot(XC(3:6,1),YC(3:6,1),'k-','LineWidth',1)
plot(XC([1 3],1),YC([1 3],1),'k-','LineWidth',1)
plot(XC([2 6],1),YC([2 6],1),'k-','LineWidth',1)

% Sensor location coordinates and labels
plot(x,y,'ks','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');

% Damage location coordinates and labels
M1x=(x(4,1)+x(5,1))/2; M1y=(y(4,1)+y(5,1))/2;
M2x=(x(6,1)+x(7,1))/2; M2y=(y(6,1)+y(7,1))/2;
M3x=(x(11,1)+x(12,1))/2; M3y=(y(11,1)+y(12,1))/2;
plot(M1x+0.05,M1y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
plot(M2x+0.05,M2y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
plot(M3x+0.05,M3y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);

% =========================================================================
plot(x,y,'rd','MarkerIndices',[bestIndicesSSCADPRFSGA2{NumSen,1}],'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7);  
plot(x,y,'rd','MarkerIndices',bestIndicesSSCADPRFSGA2{NumSen,1}(1,DelNumSSCADPRGAFS(1,NumSen)),'MarkerEdgeColor',colors(1,:),'MarkerFaceColor',colors(1,:),'MarkerSize',7); 

set(gca,'XTick',[], 'YTick', [])
axis equal
axis([-0.5 8 -0.5 1.5])
hold off

%% ============================= SSC-ADPR-FSR =============================

sf = 180;
gap = [0 0]; % [vertical,horizontal]
marg_h = [0.06 0.05]; % [lower uppper]
marg_w = [0.02 0.02]; % [left right]
subplot = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
figure('rend','painters','pos',[100 100 sf*3.4 sf])

% the outline of the wing 
figure(12)
subplot(1, 1, 1)
hold on
XC(1,1) = x(1,1)-0.15; XC(2,1) = x(12,1)+0.15;
YC(1,1) = y(1,1)-0.25; YC(2,1) = y(12,1)-0.25;
plot(XC([1 2],1),YC([1 2],1),'k-','LineWidth',1);

XC(3,1) = x(25,1)-0.15; XC(4,1) = x(31,1); XC(5,1) = x(32,1); XC(6,1) = x(36,1)+0.15;
YC(3,1) = y(25,1)+0.25; YC(4,1) = y(31,1)+0.25; YC(5,1) = y(32,1)+0.25; YC(6,1) = y(36,1)+0.25;

plot(XC(3:6,1),YC(3:6,1),'k-','LineWidth',1)
plot(XC([1 3],1),YC([1 3],1),'k-','LineWidth',1)
plot(XC([2 6],1),YC([2 6],1),'k-','LineWidth',1)

% Sensor location coordinates and labels
plot(x,y,'ks','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k');

% Damage location coordinates and labels
M1x=(x(4,1)+x(5,1))/2; M1y=(y(4,1)+y(5,1))/2;
M2x=(x(6,1)+x(7,1))/2; M2y=(y(6,1)+y(7,1))/2;
M3x=(x(11,1)+x(12,1))/2; M3y=(y(11,1)+y(12,1))/2;
plot(M1x+0.05,M1y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
plot(M2x+0.05,M2y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);
plot(M3x+0.05,M3y,'go','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);

% =========================================================================
plot(x,y,'rd','MarkerIndices',[bestIndicesSSCADPRFSRGA2{NumSen2,1}],'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7);  
plot(x,y,'rd','MarkerIndices',bestIndicesSSCADPRFSRGA2{NumSen2,1}(1,DelNumSSCADPRGAFSR(1,NumSen2)),'MarkerEdgeColor',colors(5,:),'MarkerFaceColor',colors(5,:),'MarkerSize',7); 
plot(x+0.1,y,'rd','MarkerIndices',bestIndicesSSCADPRFSRGA2{NumSen2,1}(1,RepNumSSCADPRGAFSR(1,NumSen2)),'MarkerEdgeColor',colors(7,:),'MarkerFaceColor',colors(7,:),'MarkerSize',7);  

set(gca,'XTick',[], 'YTick', [])
axis equal
axis([-0.5 8 -0.5 1.5])
hold off

