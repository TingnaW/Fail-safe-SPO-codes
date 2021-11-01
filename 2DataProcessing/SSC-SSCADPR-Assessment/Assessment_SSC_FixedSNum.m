clear;
colors = lines(5);
FS = 15;
FS2 = 14;

%% =========================== Load results ===============================
% addpath('add the location of the following files')  
addpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/2DataProcessing/SSC')

load OptResSSCExh_Mode2.mat bestIndicesSSCEH SSCEH optimalSSCEH
load OptResSSCExhFS2.mat bestIndicesSSCFSEH2 optimalSSCFSEH2
load OptResSSCExhFSR2.mat bestIndicesSSCFSREH2 optimalSSCFSREH2

load OptResSSCGA.mat bestIndicesSSCGA optimalSSCGA
load OptResSSCGAFS2.mat bestIndicesSSCFSGA2 optimalSSCFSGA2
load OptResSSCGAFSR2.mat bestIndicesSSCFSRGA2 optimalSSCFSRGA2

load Coe4ModesNorDam.mat X Y

MNum = 2; % Select the mode shape, 1:{1,1}; 2:{1,2}; 3{1,3}; 4{1,4};
CoeModeS = 10^6*X{1,MNum}; 

% rmpath('add the location of these files')
rmpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/2DataProcessing/SSC')

%% *************************** Two Sensors ********************************
NumSen=2;
NumSen2 = NumSen-1;

figure(13)
histogram(SSCEH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Exh
XLim =[optimalSSCEH(NumSen,1) optimalSSCEH(NumSen,1)];YLim =[0 7*10];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[optimalSSCFSEH2(NumSen,1) optimalSSCFSEH2(NumSen,1)];YLim =[0 7*10];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% GA
XLim =optimalSSCGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalSSCFSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
% legend([p1,p2,p4,p5],'Exhaustive search (ES)','Improved fail-safe ES','Genetic algorithm (GA)','Improved fail-safe GA','FontSize',FS2,'FontName','times')

%% *********************** After one sensor fail + Worse Case *************
% Exh
C =[];C = bestIndicesSSCEH{NumSen,1};
[FSSSCExh(1,NumSen), DelNumExh(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% ExhFS
C =[];C = bestIndicesSSCFSEH2{NumSen,1};
[FSSSCExhFS(1,NumSen), DelNumExhFS(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% GA
C =[];C = bestIndicesSSCGA{NumSen,1};
[FSSSCGA(1,NumSen), DelNumGA(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% GAFS 
C =[];C = bestIndicesSSCFSGA2{NumSen,1};
[FSSSCGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

figure(14) 
histogram(SSCEH{1,NumSen2},30,'FaceColor',colors(5,:)) 
hold on 

% Exh
XLim =[FSSSCExh(1,NumSen) FSSSCExh(1,NumSen)];YLim = [0 5];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[FSSSCExhFS(1,NumSen) FSSSCExhFS(1,NumSen)];YLim = [0 5];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% GA
XLim =FSSSCGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSSSCGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Three Sensors ******************************
NumSen=3;
NumSen2 = NumSen-1;

figure(1)
histogram(SSCEH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Exh
XLim =[optimalSSCEH(NumSen,1) optimalSSCEH(NumSen,1)];YLim =[0 6*10^2];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[optimalSSCFSEH2(NumSen,1) optimalSSCFSEH2(NumSen,1)];YLim =[0 6*10^2];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[optimalSSCFSREH2(NumSen2,1) optimalSSCFSREH2(NumSen2,1)];YLim =[0 6*10^2];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalSSCGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalSSCFSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =optimalSSCFSRGA2(NumSen2,1);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
% legend([p1,p2,p3,p4,p5,p6],'Exhaustive search (ES)','Improved fail-safe ES','Improved fail-safe with redundancy ES',...
%                  'Genetic algorithm (GA)','Improved fail-safe GA','Improved fail-safe with redundancy GA','FontSize',FS2,'FontName','times')

%% *********************** After one sensor fail + Worse Case *************
% Exh
C =[];C = bestIndicesSSCEH{NumSen,1};
[FSSSCExh(1,NumSen), DelNumExh(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% ExhFS
C =[];C = bestIndicesSSCFSEH2{NumSen,1};
[FSSSCExhFS(1,NumSen), DelNumExhFS(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% ExhFSR
C =[];C = bestIndicesSSCFSREH2{NumSen2,1};
[FSSSCExhFSR(1,NumSen2), DelNumExhFSR(1,NumSen2),RepNumExhFSR(1,NumSen2)]= SSC_FucFSR(NumSen2,C,CoeModeS,Y);

% GA
C =[];C = bestIndicesSSCGA{NumSen,1};
[FSSSCGA(1,NumSen), DelNumGA(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% GAFS 
C =[];C = bestIndicesSSCFSGA2{NumSen,1};
[FSSSCGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSRGA2{NumSen2,1};
[FSSSCGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= SSC_FucFSR(NumSen2,C,CoeModeS,Y);

figure(2) 
histogram(SSCEH{1,NumSen-1},30,'FaceColor',colors(5,:)) 
hold on 

% Exh
XLim =[FSSSCExh(1,NumSen) FSSSCExh(1,NumSen)];YLim = [0 6*10^1];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[FSSSCExhFS(1,NumSen) FSSSCExhFS(1,NumSen)];YLim = [0 6*10^1];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[FSSSCExhFSR(1,NumSen2) FSSSCExhFSR(1,NumSen2)];YLim = [0 6*10^1];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSSSCGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSSSCGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =FSSSCGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Four Sensors *******************************
NumSen=4;
NumSen2 = NumSen-1;

figure(3)
histogram(SSCEH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Exh
XLim =[optimalSSCEH(NumSen,1) optimalSSCEH(NumSen,1)];YLim =[0 7*10^3];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[optimalSSCFSEH2(NumSen,1) optimalSSCFSEH2(NumSen,1)];YLim =[0 7*10^3];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[optimalSSCFSREH2(NumSen2,1) optimalSSCFSREH2(NumSen2,1)];YLim =[0 7*10^3];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalSSCGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalSSCFSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =optimalSSCFSRGA2(NumSen2,1);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
% legend([p1,p2,p3],'Exhaustive search (ES)','Improved fail-safe ES','Improved fail-safe with redundancy ES','FontSize',FS2,'FontName','times')

%% *********************** After one sensor fail + Worse Case *************
% Exh
C =[];C = bestIndicesSSCEH{NumSen,1};
[FSSSCExh(1,NumSen), DelNumExh(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% ExhFS
C =[];C = bestIndicesSSCFSEH2{NumSen,1};
[FSSSCExhFS(1,NumSen), DelNumExhFS(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% ExhFSR
C =[];C = bestIndicesSSCFSREH2{NumSen2,1};
[FSSSCExhFSR(1,NumSen2), DelNumExhFSR(1,NumSen2),RepNumExhFSR(1,NumSen2)]= SSC_FucFSR(NumSen2,C,CoeModeS,Y);

% GA
C =[];C = bestIndicesSSCGA{NumSen,1};
[FSSSCGA(1,NumSen), DelNumGA(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% GAFS 
C =[];C = bestIndicesSSCFSGA2{NumSen,1};
[FSSSCGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSRGA2{NumSen2,1};
[FSSSCGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= SSC_FucFSR(NumSen2,C,CoeModeS,Y);

figure(4) 
histogram(SSCEH{1,NumSen-1},30,'FaceColor',colors(5,:)) 
hold on 

% Exh
XLim =[FSSSCExh(1,NumSen) FSSSCExh(1,NumSen)];YLim = [0 8*10^2];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[FSSSCExhFS(1,NumSen) FSSSCExhFS(1,NumSen)];YLim = [0 8*10^2];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[FSSSCExhFSR(1,NumSen2) FSSSCExhFSR(1,NumSen2)];YLim = [0 8*10^2];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSSSCGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSSSCGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =FSSSCGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
% legend([p4,p5,p6],'Genetic algorithm (GA)','Improved fail-safe GA','Improved fail-safe with redundancy GA','FontSize',FS2,'FontName','times')

%% *************************** Five Sensors *******************************
NumSen=5;
NumSen2 = NumSen-1;

figure(5)
histogram(SSCEH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Exh
XLim =[optimalSSCEH(NumSen,1) optimalSSCEH(NumSen,1)];YLim =[0 5*10^4];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[optimalSSCFSEH2(NumSen,1) optimalSSCFSEH2(NumSen,1)];YLim =[0 5*10^4];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[optimalSSCFSREH2(NumSen2,1) optimalSSCFSREH2(NumSen2,1)];YLim =[0 5*10^4];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalSSCGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalSSCFSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =optimalSSCFSRGA2(NumSen2,1);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case *************
% Exh
C =[];C = bestIndicesSSCEH{NumSen,1};
[FSSSCExh(1,NumSen), DelNumExh(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% ExhFS
C =[];C = bestIndicesSSCFSEH2{NumSen,1};
[FSSSCExhFS(1,NumSen), DelNumExhFS(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% ExhFSR
C =[];C = bestIndicesSSCFSREH2{NumSen2,1};
[FSSSCExhFSR(1,NumSen2), DelNumExhFSR(1,NumSen2),RepNumExhFSR(1,NumSen2)]= SSC_FucFSR(NumSen2,C,CoeModeS,Y);

% GA
C =[];C = bestIndicesSSCGA{NumSen,1};
[FSSSCGA(1,NumSen), DelNumGA(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% GAFS 
C =[];C = bestIndicesSSCFSGA2{NumSen,1};
[FSSSCGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSRGA2{NumSen2,1};
[FSSSCGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= SSC_FucFSR(NumSen2,C,CoeModeS,Y);

figure(6) 
histogram(SSCEH{1,NumSen-1},30,'FaceColor',colors(5,:)) 
hold on 

% Exh
XLim =[FSSSCExh(1,NumSen) FSSSCExh(1,NumSen)];YLim = [0 7*10^3];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[FSSSCExhFS(1,NumSen) FSSSCExhFS(1,NumSen)];YLim = [0 7*10^3];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[FSSSCExhFSR(1,NumSen2) FSSSCExhFSR(1,NumSen2)];YLim = [0 7*10^3];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSSSCGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSSSCGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =FSSSCGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Six Sensors ********************************
NumSen=6;
NumSen2 = NumSen-1;

figure(7)
histogram(SSCEH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Exh
XLim =[optimalSSCEH(NumSen,1) optimalSSCEH(NumSen,1)];YLim =[0 3*10^5];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[optimalSSCFSEH2(NumSen,1) optimalSSCFSEH2(NumSen,1)];YLim =[0 3*10^5];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[optimalSSCFSREH2(NumSen2,1) optimalSSCFSREH2(NumSen2,1)];YLim =[0 3*10^5];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalSSCGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalSSCFSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =optimalSSCFSRGA2(NumSen2,1);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case *************
% Exh
C =[];C = bestIndicesSSCEH{NumSen,1};
[FSSSCExh(1,NumSen), DelNumExh(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% ExhFS
C =[];C = bestIndicesSSCFSEH2{NumSen,1};
[FSSSCExhFS(1,NumSen), DelNumExhFS(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% ExhFSR
C =[];C = bestIndicesSSCFSREH2{NumSen2,1};
[FSSSCExhFSR(1,NumSen2), DelNumExhFSR(1,NumSen2),RepNumExhFSR(1,NumSen2)]= SSC_FucFSR(NumSen2,C,CoeModeS,Y);

% GA
C =[];C = bestIndicesSSCGA{NumSen,1};
[FSSSCGA(1,NumSen), DelNumGA(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% GAFS 
C =[];C = bestIndicesSSCFSGA2{NumSen,1};
[FSSSCGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSRGA2{NumSen2,1};
[FSSSCGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= SSC_FucFSR(NumSen2,C,CoeModeS,Y);

figure(8) 
histogram(SSCEH{1,NumSen-1},30,'FaceColor',colors(5,:)) 
hold on 

% Exh
XLim =[FSSSCExh(1,NumSen) FSSSCExh(1,NumSen)];YLim = [0 5*10^4];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[FSSSCExhFS(1,NumSen) FSSSCExhFS(1,NumSen)];YLim = [0 5*10^4];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[FSSSCExhFSR(1,NumSen2) FSSSCExhFSR(1,NumSen2)];YLim = [0 5*10^4];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSSSCGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSSSCGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =FSSSCGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Seven Sensors ******************************
NumSen=7;
NumSen2 = NumSen-1;

figure(9)
hold on 

% ExhFSR
XLim =[optimalSSCFSREH2(NumSen2,1) optimalSSCFSREH2(NumSen2,1)];YLim =[0 2*10^5];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalSSCGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalSSCFSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =optimalSSCFSRGA2(NumSen2,1);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case *************
% ExhFSR
C =[];C = bestIndicesSSCFSREH2{NumSen2,1};
[FSSSCExhFSR(1,NumSen2), DelNumExhFSR(1,NumSen2),RepNumExhFSR(1,NumSen2)]= SSC_FucFSR(NumSen2,C,CoeModeS,Y);

% GA
C =[];C = bestIndicesSSCGA{NumSen,1};
[FSSSCGA(1,NumSen), DelNumGA(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% GAFS 
C =[];C = bestIndicesSSCFSGA2{NumSen,1};
[FSSSCGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSRGA2{NumSen2,1};
[FSSSCGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= SSC_FucFSR(NumSen2,C,CoeModeS,Y);

figure(10) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
hold on 

% ExhFSR
XLim =[FSSSCExhFSR(1,NumSen2) FSSSCExhFSR(1,NumSen2)];YLim = [0 2*10^5];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSSSCGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSSSCGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =FSSSCGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Eight Sensor locations *********************
NumSen=8;
NumSen2 = NumSen-1;

figure(11)
hold on 
% GA
XLim =[optimalSSCGA(NumSen,1) optimalSSCGA(NumSen,1)];YLim = [0 2*10^5];
p2=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% GAR: Fail-safe
XLim =[optimalSSCFSGA2(NumSen,1) optimalSSCFSGA2(NumSen,1)];YLim = [0 2*10^5];
p3=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',12)

% GAR: Fail-safe with redundancy
XLim =[optimalSSCFSRGA2(NumSen2,1) optimalSSCFSRGA2(NumSen2,1)];YLim = [0 2*10^5];
p4=plot(XLim,YLim,'--','Color',colors(3,:),'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',FS)

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case *************
% GA
C =[];C = bestIndicesSSCGA{NumSen,1};
[FSSSCGA(1,NumSen), DelNumGA(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% GAR 
C =[];C = bestIndicesSSCFSGA2{NumSen,1};
[FSSSCGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= SSC_FucFS(NumSen,C,CoeModeS,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSRGA2{NumSen2,1};
[FSSSCGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= SSC_FucFSR(NumSen2,C,CoeModeS,Y);

figure(12) 
hold on 

% GA
XLim =[FSSSCGA(1,NumSen) FSSSCGA(1,NumSen)];YLim = [0 2*10^5];
p2=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% GAR: Fail-safe
XLim =[FSSSCGAFS(1,NumSen) FSSSCGAFS(1,NumSen)];YLim = [0 2*10^5];
p3=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',12)

% GAR: Fail-safe with redundancy
XLim =[FSSSCGAFSR(1,NumSen2) FSSSCGAFSR(1,NumSen2)];YLim = [0 2*10^5];
p4=plot(XLim,YLim,'--','Color',colors(3,:),'LineWidth',2);
xlabel('SSC');ylabel('Frequency');set(gca,'FontSize',FS)

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Nine Sensor locations ********************** 
NumSen=9;
NumSen2 = NumSen-1;
% GAFSR 
C =[];C = bestIndicesSSCFSRGA2{NumSen2,1};
[FSSSCGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= SSC_FucFSR(NumSen2,C,CoeModeS,Y);

%% ======== Save the results  =============================================
save SSCGASensorDistribution.mat bestIndicesSSCGA DelNumGA bestIndicesSSCFSGA2 DelNumGAFS bestIndicesSSCFSRGA2 DelNumGAFSR RepNumGAFSR

%% ======== Function 1  ===================================================
function [SSCFS, I]= SSC_FucFS(NumSen,C,CoeModeS,Y)
CoeSelM =[]; 
    for i=1:NumSen
        CoeSelM = [CoeSelM CoeModeS(:,C(1,i))];  
    end

    for i=1:NumSen
        CoeSelMMid = CoeSelM;
        CoeSelMMid(:,i) = [];
        [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
        MEdMid(1,i)= sum(CCC.^2);   
    end    
[SSCFS, I]= min(MEdMid);
end

%% ======== Function 2  ===================================================
function [SSCFSR, Idel, IRep]= SSC_FucFSR(NumSen,C,CoeModeS,Y)
CoeSecM =[]; 
    for i=1:NumSen
        CoeSecM = [CoeSecM CoeModeS(:,C(1,i))];  
    end

    for i=1:NumSen
        CoeSelMMid = CoeSecM;
        CoeSelMMid(:,i) = [];
        [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
        MEdMid(1,i)= sum(CCC.^2); 
    end    
[~,IRep] = min(MEdMid);
MEdMidIni = MEdMid;
MEdMid(IRep) = [];   
SSCFSR = min(MEdMid);
Idel = find(MEdMidIni==SSCFSR,1, 'first');
end
