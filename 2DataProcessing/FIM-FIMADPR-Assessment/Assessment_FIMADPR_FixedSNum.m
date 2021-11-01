clear;
colors = lines(5);
FS = 15;
FS2 = 14;

%% =========================== Load results ===============================
% addpath('add the location of the following files')  
addpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/2DataProcessing/FIMADPR')

load OptResFIMADPRExh.mat bestIndicesFIMADPREH FIMADPREH optimalFIMADPREH
load OptResFIMADPRExhFS2.mat bestIndicesFIMADPR_FSEH2 optimalFIMADPR_FSEH2
load OptResFIMADPRExhFSR2.mat bestIndicesFIMADPR_FSREH2 optimalFIMADPR_FSREH2

load OptResFIMADPRGA.mat bestIndicesFIMADPRGA optimalFIMADPRGA
load OptResFIMADPRGAFS2.mat bestIndicesFIMADPR_FSGA2 optimalFIMADPR_FSGA2
load OptResFIMADPRGAFSR2.mat bestIndicesFIMADPR_FSRGA2 optimalFIMADPR_FSRGA2

load Coe4Modes.mat Coe4Modes NF
Coe4ModeM =10^3*Coe4Modes{1,4};

% rmpath('add the location of these files')
rmpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/2DataProcessing/FIMADPR')

%% *************************** Four Sensors *******************************
NumSen=4;
NumSen2 = NumSen-1;

figure(9)
histogram(FIMADPREH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Ehx
XLim =[optimalFIMADPREH(NumSen,1) optimalFIMADPREH(NumSen,1)];YLim = [0 5*10^4];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[optimalFIMADPR_FSEH2(NumSen,1) optimalFIMADPR_FSEH2(NumSen,1)];YLim = [0 5*10^4];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% GA
XLim =optimalFIMADPRGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFIMADPR_FSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
% legend([p1,p2,p4,p5],'Exhaustive search (ES)','Improved fail-safe ES','Genetic algorithm (GA)','Improved fail-safe GA','FontSize',FS2,'FontName','times')

%% *********************** After one sensor fail + Worse Case *************
% Exh
C =[];C = bestIndicesFIMADPREH{NumSen,1};
[FSFIMADPRExh(1,NumSen), DelNumADPRExh(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% Exh
C =[];C = bestIndicesFIMADPR_FSEH2{NumSen,1};
[FSFIMADPRExhFS(1,NumSen), DelNumADPRExhFS(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% GA
C =[];C = bestIndicesFIMADPRGA{NumSen,1};
[FSFIMADPRGA(1,NumSen), DelNumADPRGA(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% GAFS 
C =[];C = bestIndicesFIMADPR_FSGA2{NumSen,1};
[FSFIMADPRGAFS(1,NumSen), DelNumADPRGAFS(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

figure(10) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
histogram(FIMADPREH{1,NumSen2},30,'FaceColor',colors(5,:)) 
hold on 

% Ehx
XLim =[FSFIMADPRExh(1,NumSen) FSFIMADPRExh(1,NumSen)];YLim = [0 7*10^3];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMADPRExhFS(1,NumSen) FSFIMADPRExhFS(1,NumSen)];YLim = [0 7*10^3];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% GA
XLim =FSFIMADPRGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSFIMADPRGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Five Sensors *******************************
NumSen=5;
NumSen2 = NumSen-1;

figure(1)
histogram(FIMADPREH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Ehx
XLim =[optimalFIMADPREH(NumSen,1) optimalFIMADPREH(NumSen,1)];YLim = [0 3*10^5];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[optimalFIMADPR_FSEH2(NumSen,1) optimalFIMADPR_FSEH2(NumSen,1)];YLim = [0 3*10^5];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[optimalFIMADPR_FSREH2(NumSen2,1) optimalFIMADPR_FSREH2(NumSen2,1)];YLim = [0 3*10^5];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalFIMADPRGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFIMADPR_FSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =optimalFIMADPR_FSRGA2(NumSen2,1);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
% legend([p3,p6],'Improved fail-safe with redundancy ES','Improved fail-safe with redundancy GA','FontSize',FS2,'FontName','times')

%% *********************** After one sensor fail + Worse Case *************
% Exh
C =[];C = bestIndicesFIMADPREH{NumSen,1};
[FSFIMADPRExh(1,NumSen), DelNumADPRExh(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% Exh
C =[];C = bestIndicesFIMADPR_FSEH2{NumSen,1};
[FSFIMADPRExhFS(1,NumSen), DelNumADPRExhFS(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% ExhFSR
C =[];C = bestIndicesFIMADPR_FSREH2{NumSen2,1};
[FSFIMADPRExhFSR(1,NumSen2), DelNumADPRExhFSR(1,NumSen2),RepNumADPRExhFSR(1,NumSen2)]= FIMADPR_FucFSR(NumSen2,C,Coe4ModeM,NF);

% GA
C =[];C = bestIndicesFIMADPRGA{NumSen,1};
[FSFIMADPRGA(1,NumSen), DelNumADPRGA(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% GAFS 
C =[];C = bestIndicesFIMADPR_FSGA2{NumSen,1};
[FSFIMADPRGAFS(1,NumSen), DelNumADPRGAFS(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% GAFSR
C =[];C = bestIndicesFIMADPR_FSRGA2{NumSen2,1};
[FSFIMADPRGAFSR(1,NumSen2), DelNumADPRGAFSR(1,NumSen2),RepNumADPRGAFSR(1,NumSen2)]= FIMADPR_FucFSR(NumSen2,C,Coe4ModeM,NF);

figure(2) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
histogram(FIMADPREH{1,NumSen2},30,'FaceColor',colors(5,:)) % ^^^^^^^^^^^^^^
hold on 

% Ehx
XLim =[FSFIMADPRExh(1,NumSen) FSFIMADPRExh(1,NumSen)];YLim = [0 5*10^4];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMADPRExhFS(1,NumSen) FSFIMADPRExhFS(1,NumSen)];YLim = [0 5*10^4];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMADPRExhFSR(1,NumSen2) FSFIMADPRExhFSR(1,NumSen2)];YLim = [0 5*10^4];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSFIMADPRGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSFIMADPRGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =FSFIMADPRGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Six Sensors ********************************
NumSen=6;
NumSen2 = NumSen-1;

figure(3)
histogram(FIMADPREH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Ehx
XLim =[optimalFIMADPREH(NumSen,1) optimalFIMADPREH(NumSen,1)];YLim = [0 12*10^5];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[optimalFIMADPR_FSEH2(NumSen,1) optimalFIMADPR_FSEH2(NumSen,1)];YLim = [0 12*10^5];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[optimalFIMADPR_FSREH2(NumSen2,1) optimalFIMADPR_FSREH2(NumSen2,1)];YLim = [0 12*10^5];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalFIMADPRGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFIMADPR_FSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =optimalFIMADPR_FSRGA2(NumSen2,1);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case *************
% Exh
C =[];C = bestIndicesFIMADPREH{NumSen,1};
[FSFIMADPRExh(1,NumSen), DelNumADPRExh(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% Exh
C =[];C = bestIndicesFIMADPR_FSEH2{NumSen,1};
[FSFIMADPRExhFS(1,NumSen), DelNumADPRExhFS(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% ExhFSR
C =[];C = bestIndicesFIMADPR_FSREH2{NumSen2,1};
[FSFIMADPRExhFSR(1,NumSen2), DelNumADPRExhFSR(1,NumSen2),RepNumADPRExhFSR(1,NumSen2)]= FIMADPR_FucFSR(NumSen2,C,Coe4ModeM,NF);

% GA
C =[];C = bestIndicesFIMADPRGA{NumSen,1};
[FSFIMADPRGA(1,NumSen), DelNumADPRGA(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% GAFS 
C =[];C = bestIndicesFIMADPR_FSGA2{NumSen,1};
[FSFIMADPRGAFS(1,NumSen), DelNumADPRGAFS(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% GAFSR
C =[];C = bestIndicesFIMADPR_FSRGA2{NumSen2,1};
[FSFIMADPRGAFSR(1,NumSen2), DelNumADPRGAFSR(1,NumSen2),RepNumADPRGAFSR(1,NumSen2)]= FIMADPR_FucFSR(NumSen2,C,Coe4ModeM,NF);

figure(4) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
histogram(FIMADPREH{1,NumSen2},30,'FaceColor',colors(5,:)) 
hold on 

% Ehx
XLim =[FSFIMADPRExh(1,NumSen) FSFIMADPRExh(1,NumSen)];YLim = [0 3*10^5];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMADPRExhFS(1,NumSen) FSFIMADPRExhFS(1,NumSen)];YLim = [0 3*10^5];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMADPRExhFSR(1,NumSen2) FSFIMADPRExhFSR(1,NumSen2)];YLim = [0 3*10^5];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSFIMADPRGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSFIMADPRGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =FSFIMADPRGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Seven Sensors ******************************
NumSen=7;
NumSen2 = NumSen-1;

figure(5)
hold on 

% ExhFSR
XLim =[optimalFIMADPR_FSREH2(NumSen2,1) optimalFIMADPR_FSREH2(NumSen2,1)];YLim = [0 3*10^5];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalFIMADPRGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFIMADPR_FSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =optimalFIMADPR_FSRGA2(NumSen2,1);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case *************
% ExhFSR
C =[];C = bestIndicesFIMADPR_FSREH2{NumSen2,1};
[FSFIMADPRExhFSR(1,NumSen2), DelNumADPRExhFSR(1,NumSen2),RepNumADPRExhFSR(1,NumSen2)]= FIMADPR_FucFSR(NumSen2,C,Coe4ModeM,NF);

% GA
C =[];C = bestIndicesFIMADPRGA{NumSen,1};
[FSFIMADPRGA(1,NumSen), DelNumADPRGA(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% GAFS 
C =[];C = bestIndicesFIMADPR_FSGA2{NumSen,1};
[FSFIMADPRGAFS(1,NumSen), DelNumADPRGAFS(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% GAFSR
C =[];C = bestIndicesFIMADPR_FSRGA2{NumSen2,1};
[FSFIMADPRGAFSR(1,NumSen2), DelNumADPRGAFSR(1,NumSen2),RepNumADPRGAFSR(1,NumSen2)]= FIMADPR_FucFSR(NumSen2,C,Coe4ModeM,NF);

figure(6) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
hold on 

% EhxFS
XLim =[FSFIMADPRExhFSR(1,NumSen2) FSFIMADPRExhFSR(1,NumSen2)];YLim = [0 3*10^5];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSFIMADPRGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSFIMADPRGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =FSFIMADPRGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Eight Sensors ******************************
NumSen=8;
NumSen2 = NumSen-1;

figure(7)
hold on 

% GA
XLim =[optimalFIMADPRGA(NumSen,1) optimalFIMADPRGA(NumSen,1)];YLim = [0 3*10^5];
p2=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% GAFS
XLim =[optimalFIMADPR_FSGA2(NumSen,1) optimalFIMADPR_FSGA2(NumSen,1)];YLim = [0 3*10^5];
p3=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);
xlabel('DFIM weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS)

% GAFSR
XLim =[optimalFIMADPR_FSRGA2(NumSen2,1) optimalFIMADPR_FSRGA2(NumSen2,1)];YLim = [0 3*10^5];
p4=plot(XLim,YLim,'--','Color',colors(3,:),'LineWidth',2);
xlabel('DFIM weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS)

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p2,p3,p4],'Genetic algorithm','Improved Fail-safe (GA)','Improved Fail-safe with redundancy (GA)','FontSize',FS2)

%% *********************** After one sensor fail + Worse Case *************

% GA
C =[];C = bestIndicesFIMADPRGA{NumSen,1};
[FSFIMADPRGA(1,NumSen), DelNumADPRGA(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% GAFS 
C =[];C = bestIndicesFIMADPR_FSGA2{NumSen,1};
[FSFIMADPRGAFS(1,NumSen), DelNumADPRGAFS(1,NumSen)]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF);

% GAFSR
C =[];C = bestIndicesFIMADPR_FSRGA2{NumSen2,1};
[FSFIMADPRGAFSR(1,NumSen2), DelNumADPRGAFSR(1,NumSen2),RepNumADPRGAFSR(1,NumSen2)]= FIMADPR_FucFSR(NumSen2,C,Coe4ModeM,NF);

figure(8) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
hold on 

% GA
XLim =[FSFIMADPRGA(1,NumSen) FSFIMADPRGA(1,NumSen)];YLim = [0 3*10^5];
p2=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% GAFS
XLim =[FSFIMADPRGAFS(1,NumSen) FSFIMADPRGAFS(1,NumSen)];YLim = [0 3*10^5];
p3=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);
xlabel('DFIM weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS)

% GAFSR
XLim =[FSFIMADPRGAFSR(1,NumSen2) FSFIMADPRGAFSR(1,NumSen2)];YLim = [0 3*10^5];
p4=plot(XLim,YLim,'--','Color',colors(3,:),'LineWidth',2);
xlabel('DFIM weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS)

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p2,p3,p4],'Genetic algorithm','Improved Fail-safe (GA)','Improved Fail-safe with redundancy (GA)','FontSize',FS2)

%% *************************** Nine Sensor locations **********************
% GAFSR
NumSen=9;
NumSen2 = NumSen-1;
C =[];C = bestIndicesFIMADPR_FSRGA2{NumSen2,1};
[FSFIMADPRGAFSR(1,NumSen2), DelNumADPRGAFSR(1,NumSen2),RepNumADPRGAFSR(1,NumSen2)]= FIMADPR_FucFSR(NumSen2,C,Coe4ModeM,NF);

%% ======== Save the results  =============================================
save FIMADPRGASensorDistribution.mat bestIndicesFIMADPRGA DelNumADPRGA bestIndicesFIMADPR_FSGA2 DelNumADPRGAFS bestIndicesFIMADPR_FSRGA2 DelNumADPRGAFSR RepNumADPRGAFSR

%% ======== Function 1  ===================================================
function [FIMADPRR, I]= FIMADPR_FucFS(NumSen,C,Coe4ModeM,NF)
DOFs = size(Coe4ModeM,1);
MSNum = size(NF,2);
ADPR=[];               
    for i=1:DOFs
        for k=1:MSNum
            ADPR(i,k)= Coe4ModeM(i,k)^2/(NF(1,k)*2*pi);
        end
    end
    ADPRDOF= sum(ADPR,2);
    ADPRDOF = normalize(ADPRDOF,'range');        
    
    CoeSecM =[]; ADPRM = [];
    for i=1:NumSen
        CoeSecM = [CoeSecM;Coe4ModeM(C(1,i),:)]; 
        ADPRM = [ADPRM;ADPRDOF(C(1,i),1)];
    end
    
    CoeSecMMid =[]; ADPRMMid =[];
    for i=1:NumSen
        CoeSecMMid = CoeSecM;
        CoeSecMMid(i,:) = [];
        ADPRMMid = ADPRM;
        ADPRMMid(i,:) =[];
        
        MEdMid(1,i)= det(CoeSecMMid.'*CoeSecMMid)*sum(ADPRMMid); 
    end    
[FIMADPRR, I] = min(MEdMid);
end

%% ======== Function 2  ===================================================
function [FIMADPRR, Idel, IRep]= FIMADPR_FucFSR(NumSen,C,Coe4ModeM,NF)
DOFs = size(Coe4ModeM,1);
MSNum = size(NF,2);
    for i=1:DOFs
        for k=1:MSNum
            ADPR(i,k)= Coe4ModeM(i,k)^2/(NF(1,k)*2*pi);
        end
    end
    ADPRDOF= sum(ADPR,2);
    ADPRDOF = normalize(ADPRDOF,'range');        
    
    CoeSecM =[]; ADPRM = [];
    for i=1:NumSen
        CoeSecM = [CoeSecM;Coe4ModeM(C(1,i),:)]; 
        ADPRM = [ADPRM;ADPRDOF(C(1,i),1)];
    end
    
    CoeSecMMid =[]; ADPRMMid =[];
    for i=1:NumSen
        CoeSecMMid = CoeSecM;
        CoeSecMMid(i,:) = [];
        ADPRMMid = ADPRM;
        ADPRMMid(i,:) =[];
        
        MEdMid(1,i)= det(CoeSecMMid.'*CoeSecMMid)*sum(ADPRMMid); 
    end   
[~,IRep] = min(MEdMid);
MEdMidIni = MEdMid;
MEdMid(IRep) = [];   
FIMADPRR = min(MEdMid);
Idel = find(MEdMidIni==FIMADPRR,1, 'first');
end
