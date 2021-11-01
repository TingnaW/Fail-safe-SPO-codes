clear;
colors = lines(5);
FS = 15;
FS2 = 14;

%% =========================== Load results ===============================
% addpath('add the location of the following files')  
addpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/2DataProcessing/SSCADPR')

load OptResSSCADPRExh.mat bestIndicesSSCADPREH SSCADPREH optimalSSCADPREH
load OptResSSCADPRExhFS2.mat bestIndicesSSCADPRFSEH2 optimalSSCADPRFSEH2
load OptResSSCADPRExhFSR2.mat bestIndicesSSCADPRFSREH2 optimalSSCADPRFSREH2

load OptResSSCADPRGA.mat bestIndicesSSCADPRGA optimalSSCADPRGA
load OptResSSCADPRGAFS2.mat bestIndicesSSCADPRFSGA2 optimalSSCADPRFSGA2
load OptResSSCADPRGAFSR2.mat bestIndicesSSCADPRFSRGA2 optimalSSCADPRFSRGA2

load Coe4ModesNorDam.mat X Y NF

MNum = 2; % Select the mode shape, 1:{1,1}; 2:{1,2}; 3{1,3}; 4{1,4};
CoeModeS = 10^6*X{1,MNum}; 

% rmpath('add the location of these files')
rmpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/2DataProcessing/SSCADPR')

%% *************************** Two Sensors ********************************
NumSen=2;
NumSen2 = NumSen-1;

figure(13)
histogram(SSCADPREH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Exh
XLim =[optimalSSCADPREH(NumSen,1) optimalSSCADPREH(NumSen,1)];YLim = [0 2*10^2];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[optimalSSCADPRFSEH2(NumSen,1) optimalSSCADPRFSEH2(NumSen,1)];YLim = [0 2*10^2];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% GA
XLim =optimalSSCADPRGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS: Fail-safe
XLim =optimalSSCADPRFSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p1,p2,p4,p5],'Exhaustive search (ES)','Improved fail-safe ES','Genetic algorithm (GA)','Improved fail-safe GA','FontSize',FS2,'FontName','times')

%% *********************** After one sensor fail + Worse Case *************
% Exh
C =[];C = bestIndicesSSCADPREH{NumSen,1};
[FSSSCADPRExh(1,NumSen), DelNumExh(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% ExhFS
C =[];C = bestIndicesSSCADPRFSEH2{NumSen,1};
[FSSSCADPRExhFS(1,NumSen), DelNumExhFS(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% GA
C =[];C = bestIndicesSSCADPRGA{NumSen,1};
[FSSSCADPRGA(1,NumSen), DelNumGA(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% GAFS 
C =[];C = bestIndicesSSCADPRFSGA2{NumSen,1};
[FSSSCADPRGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

figure(14) 
histogram(SSCADPREH{1,NumSen2},30,'FaceColor',colors(5,:)) 
hold on 

% Exh
XLim =[FSSSCADPRExh(1,NumSen) FSSSCADPRExh(1,NumSen)];YLim = [0 2*10^1];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[FSSSCADPRExhFS(1,NumSen) FSSSCADPRExhFS(1,NumSen)];YLim = [0 2*10^1];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% GA
XLim =FSSSCADPRGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',10,'LineWidth',2);

% GAR: Fail-safe
XLim =FSSSCADPRGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',10,'LineWidth',2);
xlabel('SSC weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p1,p2,p4,p5],'Exhaustive search (ES)','Improved fail-safe ES','Genetic algorithm (GA)','Improved fail-safe GA','FontSize',FS2,'FontName','times')

%% *************************** Three Sensors ******************************
NumSen=3;
NumSen2 = NumSen-1;

figure(1)
histogram(SSCADPREH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Exh
XLim =[optimalSSCADPREH(NumSen,1) optimalSSCADPREH(NumSen,1)];YLim = [0 2*10^3];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[optimalSSCADPRFSEH2(NumSen,1) optimalSSCADPRFSEH2(NumSen,1)];YLim = [0 2*10^3];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[optimalSSCADPRFSREH2(NumSen2,1) optimalSSCADPRFSREH2(NumSen2,1)];YLim = [0 2*10^3];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalSSCADPRGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS: Fail-safe
XLim =optimalSSCADPRFSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =optimalSSCADPRFSRGA2(NumSen2,1);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p1,p2,p3,p4,p5,p6],'Exhaustive search (ES)','Improved fail-safe ES','Improved fail-safe with redundancy ES',...
                 'Genetic algorithm (GA)','Improved fail-safe GA','Improved fail-safe with redundancy GA','FontSize',FS2,'FontName','times')

%% *********************** After one sensor fail + Worse Case *************
% Exh
C =[];C = bestIndicesSSCADPREH{NumSen,1};
[FSSSCADPRExh(1,NumSen), DelNumExh(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% ExhFS
C =[];C = bestIndicesSSCADPRFSEH2{NumSen,1};
[FSSSCADPRExhFS(1,NumSen), DelNumExhFS(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% ExhFSR
C =[];C = bestIndicesSSCADPRFSREH2{NumSen2,1};
[FSSSCADPRExhFSR(1,NumSen2), DelNumExhFSR(1,NumSen2), RepNumExhFSR(1,NumSen2)]= SSCADPR_FucFSR(NumSen2,C,CoeModeS,Y,NF);

% GA
C =[];C = bestIndicesSSCADPRGA{NumSen,1};
[FSSSCADPRGA(1,NumSen), DelNumGA(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% GAR 
C =[];C = bestIndicesSSCADPRFSGA2{NumSen,1};
[FSSSCADPRGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% GAFSR 
C =[];C = bestIndicesSSCADPRFSRGA2{NumSen2,1};
[FSSSCADPRGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= SSCADPR_FucFSR(NumSen2,C,CoeModeS,Y,NF);

figure(2) 
histogram(SSCADPREH{1,NumSen-1},30,'FaceColor',colors(5,:)) 
hold on 

% Exh
XLim =[FSSSCADPRExh(1,NumSen) FSSSCADPRExh(1,NumSen)];YLim = [0 2*10^2];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[FSSSCADPRExhFS(1,NumSen) FSSSCADPRExhFS(1,NumSen)];YLim = [0 2*10^2];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[FSSSCADPRExhFSR(1,NumSen2) FSSSCADPRExhFSR(1,NumSen2)];YLim = [0 2*10^2];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSSSCADPRGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe
XLim =FSSSCADPRGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =FSSSCADPRGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Four Sensors *******************************
NumSen=4;
NumSen2 = NumSen-1;

figure(3)
histogram(SSCADPREH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Exh
XLim =[optimalSSCADPREH(NumSen,1) optimalSSCADPREH(NumSen,1)];YLim = [0 1.4*10^4];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[optimalSSCADPRFSEH2(NumSen,1) optimalSSCADPRFSEH2(NumSen,1)];YLim = [0 1.4*10^4];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[optimalSSCADPRFSREH2(NumSen2,1) optimalSSCADPRFSREH2(NumSen2,1)];YLim = [0 1.4*10^4];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalSSCADPRGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS: Fail-safe
XLim =optimalSSCADPRFSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =optimalSSCADPRFSRGA2(NumSen2,1);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
% legend([p1,p2,p3,p4,p5,p6],'Exhaustive search (ES)','Improved fail-safe ES','Improved fail-safe with redundancy ES',...
%                  'Genetic algorithm (GA)','Improved fail-safe GA','Improved fail-safe with redundancy GA','FontSize',FS2,'FontName','times')
% legend([p1,p2,p3],'Exhaustive search (ES)','Improved fail-safe ES','Improved fail-safe with redundancy ES','FontSize',FS2,'FontName','times')

%% *********************** After one sensor fail + Worse Case *************
% Exh
C =[];C = bestIndicesSSCADPREH{NumSen,1};
[FSSSCADPRExh(1,NumSen), DelNumExh(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% ExhFS
C =[];C = bestIndicesSSCADPRFSEH2{NumSen,1};
[FSSSCADPRExhFS(1,NumSen), DelNumExhFS(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% ExhFSR
C =[];C = bestIndicesSSCADPRFSREH2{NumSen2,1};
[FSSSCADPRExhFSR(1,NumSen2), DelNumExhFSR(1,NumSen2), RepNumExhFSR(1,NumSen2)]= SSCADPR_FucFSR(NumSen2,C,CoeModeS,Y,NF);

% GA
C =[];C = bestIndicesSSCADPRGA{NumSen,1};
[FSSSCADPRGA(1,NumSen), DelNumGA(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% GAR 
C =[];C = bestIndicesSSCADPRFSGA2{NumSen,1};
[FSSSCADPRGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% GAFSR 
C =[];C = bestIndicesSSCADPRFSRGA2{NumSen2,1};
[FSSSCADPRGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= SSCADPR_FucFSR(NumSen2,C,CoeModeS,Y,NF);

figure(4) 
histogram(SSCADPREH{1,NumSen-1},30,'FaceColor',colors(5,:)) 
hold on 

% Exh
XLim =[FSSSCADPRExh(1,NumSen) FSSSCADPRExh(1,NumSen)];YLim = [0 2*10^3];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[FSSSCADPRExhFS(1,NumSen) FSSSCADPRExhFS(1,NumSen)];YLim = [0 2*10^3];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[FSSSCADPRExhFSR(1,NumSen2) FSSSCADPRExhFSR(1,NumSen2)];YLim = [0 2*10^3];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSSSCADPRGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe
XLim =FSSSCADPRGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =FSSSCADPRGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
% legend([p4,p5,p6],'Genetic algorithm (GA)','Improved fail-safe GA','Improved fail-safe with redundancy GA','FontSize',FS2,'FontName','times')

%% *************************** Five Sensors *******************************
NumSen=5;
NumSen2 = NumSen-1;

figure(5)
histogram(SSCADPREH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Exh
XLim =[optimalSSCADPREH(NumSen,1) optimalSSCADPREH(NumSen,1)];YLim = [0 8*10^4];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[optimalSSCADPRFSEH2(NumSen,1) optimalSSCADPRFSEH2(NumSen,1)];YLim = [0 8*10^4];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[optimalSSCADPRFSREH2(NumSen2,1) optimalSSCADPRFSREH2(NumSen2,1)];YLim = [0 8*10^4];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalSSCADPRGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS: Fail-safe
XLim =optimalSSCADPRFSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =optimalSSCADPRFSRGA2(NumSen2,1);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case ***************
% Exh
C =[];C = bestIndicesSSCADPREH{NumSen,1};
[FSSSCADPRExh(1,NumSen), DelNumExh(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% ExhFS
C =[];C = bestIndicesSSCADPRFSEH2{NumSen,1};
[FSSSCADPRExhFS(1,NumSen), DelNumExhFS(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% ExhFSR
C =[];C = bestIndicesSSCADPRFSREH2{NumSen2,1};
[FSSSCADPRExhFSR(1,NumSen2), DelNumExhFSR(1,NumSen2), RepNumExhFSR(1,NumSen2)]= SSCADPR_FucFSR(NumSen2,C,CoeModeS,Y,NF);

% GA
C =[];C = bestIndicesSSCADPRGA{NumSen,1};
[FSSSCADPRGA(1,NumSen), DelNumGA(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% GAR 
C =[];C = bestIndicesSSCADPRFSGA2{NumSen,1};
[FSSSCADPRGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% GAFSR 
C =[];C = bestIndicesSSCADPRFSRGA2{NumSen2,1};
[FSSSCADPRGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= SSCADPR_FucFSR(NumSen2,C,CoeModeS,Y,NF);

figure(6) 
histogram(SSCADPREH{1,NumSen-1},30,'FaceColor',colors(5,:)) 
hold on 

% Exh
XLim =[FSSSCADPRExh(1,NumSen) FSSSCADPRExh(1,NumSen)];YLim = [0 1.4*10^4];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[FSSSCADPRExhFS(1,NumSen) FSSSCADPRExhFS(1,NumSen)];YLim = [0 1.4*10^4];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[FSSSCADPRExhFSR(1,NumSen2) FSSSCADPRExhFSR(1,NumSen2)];YLim = [0 1.4*10^4];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSSSCADPRGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe
XLim =FSSSCADPRGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =FSSSCADPRGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Six Sensors ********************************
NumSen=6;
NumSen2 = NumSen-1;

figure(7)
histogram(SSCADPREH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Exh
XLim =[optimalSSCADPREH(NumSen,1) optimalSSCADPREH(NumSen,1)];YLim = [0 3.5*10^5];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[optimalSSCADPRFSEH2(NumSen,1) optimalSSCADPRFSEH2(NumSen,1)];YLim = [0 3.5*10^5];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[optimalSSCADPRFSREH2(NumSen2,1) optimalSSCADPRFSREH2(NumSen2,1)];YLim = [0 3.5*10^5];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalSSCADPRGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS: Fail-safe
XLim =optimalSSCADPRFSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =optimalSSCADPRFSRGA2(NumSen2,1);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case ***************
% Exh
C =[];C = bestIndicesSSCADPREH{NumSen,1};
[FSSSCADPRExh(1,NumSen), DelNumExh(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% ExhFS
C =[];C = bestIndicesSSCADPRFSEH2{NumSen,1};
[FSSSCADPRExhFS(1,NumSen), DelNumExhFS(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% ExhFSR
C =[];C = bestIndicesSSCADPRFSREH2{NumSen2,1};
[FSSSCADPRExhFSR(1,NumSen2), DelNumExhFSR(1,NumSen2), RepNumExhFSR(1,NumSen2)]= SSCADPR_FucFSR(NumSen2,C,CoeModeS,Y,NF);

% GA
C =[];C = bestIndicesSSCADPRGA{NumSen,1};
[FSSSCADPRGA(1,NumSen), DelNumGA(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% GAR 
C =[];C = bestIndicesSSCADPRFSGA2{NumSen,1};
[FSSSCADPRGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% GAFSR 
C =[];C = bestIndicesSSCADPRFSRGA2{NumSen2,1};
[FSSSCADPRGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= SSCADPR_FucFSR(NumSen2,C,CoeModeS,Y,NF);

figure(8) 
histogram(SSCADPREH{1,NumSen-1},30,'FaceColor',colors(5,:)) 
hold on 

% Exh
XLim =[FSSSCADPRExh(1,NumSen) FSSSCADPRExh(1,NumSen)];YLim = [0 8*10^4];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% ExhFS
XLim =[FSSSCADPRExhFS(1,NumSen) FSSSCADPRExhFS(1,NumSen)];YLim = [0 8*10^4];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[FSSSCADPRExhFSR(1,NumSen2) FSSSCADPRExhFSR(1,NumSen2)];YLim = [0 8*10^4];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSSSCADPRGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe
XLim =FSSSCADPRGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =FSSSCADPRGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Seven Sensors ******************************
NumSen=7;
NumSen2 = NumSen-1;

figure(9)
hold on 

% ExhFSR
XLim =[optimalSSCADPRFSREH2(NumSen2,1) optimalSSCADPRFSREH2(NumSen2,1)];YLim = [0 2*10^2];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalSSCADPRGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS: Fail-safe
XLim =optimalSSCADPRFSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =optimalSSCADPRFSRGA2(NumSen2,1);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case *************

% ExhFSR
C =[];C = bestIndicesSSCADPRFSREH2{NumSen2,1};
[FSSSCADPRExhFSR(1,NumSen2), DelNumExhFSR(1,NumSen2), RepNumExhFSR(1,NumSen2)]= SSCADPR_FucFSR(NumSen2,C,CoeModeS,Y,NF);

% GA
C =[];C = bestIndicesSSCADPRGA{NumSen,1};
[FSSSCADPRGA(1,NumSen), DelNumGA(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% GAR 
C =[];C = bestIndicesSSCADPRFSGA2{NumSen,1};
[FSSSCADPRGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% GAFSR 
C =[];C = bestIndicesSSCADPRFSRGA2{NumSen2,1};
[FSSSCADPRGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= SSCADPR_FucFSR(NumSen2,C,CoeModeS,Y,NF);

figure(10)
hold on 

% ExhFSR
XLim =[FSSSCADPRExhFSR(1,NumSen2) FSSSCADPRExhFSR(1,NumSen2)];YLim = [0 2*10^2];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSSSCADPRGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe
XLim =FSSSCADPRGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =FSSSCADPRGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('SSC weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Eight Sensor locations *********************
NumSen=8;
NumSen2 = NumSen-1;

figure(11)
hold on 
% GA
XLim =[optimalSSCADPRGA(NumSen,1) optimalSSCADPRGA(NumSen,1)];YLim = [0 2*10^2];
p2=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% GAR: Fail-safe
XLim =[optimalSSCADPRFSGA2(NumSen,1) optimalSSCADPRFSGA2(NumSen,1)];YLim = [0 2*10^2];
p3=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =[optimalSSCADPRFSRGA2(NumSen2,1) optimalSSCADPRFSRGA2(NumSen2,1)];YLim = [0 2*10^2];
p4=plot(XLim,YLim,'--','Color',colors(3,:),'LineWidth',2);
xlabel('SSC weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS)

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case ***************

% GA
C =[];C = bestIndicesSSCADPRGA{NumSen,1};
[FSSSCADPRGA(1,NumSen), DelNumGA(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% GAR 
C =[];C = bestIndicesSSCADPRFSGA2{NumSen,1};
[FSSSCADPRGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF);

% GAFSR 
C =[];C = bestIndicesSSCADPRFSRGA2{NumSen2,1};
[FSSSCADPRGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= SSCADPR_FucFSR(NumSen2,C,CoeModeS,Y,NF);

figure(12) 
hold on 

% GA
XLim =[FSSSCADPRGA(1,NumSen) FSSSCADPRGA(1,NumSen)];YLim = [0 2*10^2];
p2=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% GAR: Fail-safe
XLim =[FSSSCADPRGAFS(1,NumSen) FSSSCADPRGAFS(1,NumSen)];YLim = [0 2*10^2];
p3=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =[FSSSCADPRGAFSR(1,NumSen2) FSSSCADPRGAFSR(1,NumSen2)];YLim = [0 2*10^2];
p4=plot(XLim,YLim,'--','Color',colors(3,:),'LineWidth',2);
xlabel('SSC weighted by ADPR');ylabel('Frequency');set(gca,'FontSize',FS)

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Nine Sensor locations ********************** 
% GAFSR 
NumSen=9;
NumSen2 = NumSen-1;
C =[];C = bestIndicesSSCADPRFSRGA2{NumSen2,1};
[FSSSCADPRGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= SSCADPR_FucFSR(NumSen2,C,CoeModeS,Y,NF);

%% ======== Save the results  =============================================
save SSCADPRGASensorDistribution.mat bestIndicesSSCADPRGA DelNumGA bestIndicesSSCADPRFSGA2 DelNumGAFS bestIndicesSSCADPRFSRGA2 DelNumGAFSR RepNumGAFSR

%% ======== Function 1  ===================================================
function [SSCFS, I]= SSCADPR_FucFS(NumSen,C,CoeModeS,Y,NF)

    NumObe = size(CoeModeS,1);
    DoFs = size(CoeModeS,2);
    ADPR=[];               % *****************
    for i=1:DoFs
        for k=1:NumObe
            ADPR(k,i)= CoeModeS(k,i)^2/(NF(k,2)*2*pi); % (NF(k,2); 2nd NF
        end
    end
    ADPRDOF= sum(ADPR,1);
    ADPRDOF = normalize(ADPRDOF,'range');  

    CoeSelM =[]; ADPRM=[];
    for i=1:NumSen
        CoeSelM = [CoeSelM CoeModeS(:,C(1,i))];  
        ADPRM = [ADPRM ADPRDOF(1,C(1,i))];
    end

    for i=1:NumSen
        CoeSelMMid = CoeSelM;
        CoeSelMMid(:,i) = [];
        ADPRMMid = ADPRM;
        ADPRMMid(:,i) = []; 
        
        [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
        MEdMid(1,i)= sum(CCC.^2)*sum(ADPRMMid);   
    end    
[SSCFS, I]= min(MEdMid);
end

%% ======== Function 2  ===================================================
function [SSCFSR, Idel, IRep]= SSCADPR_FucFSR(NumSen,C,CoeModeS,Y,NF)

    NumObe = size(CoeModeS,1);
    DoFs = size(CoeModeS,2);
    for i=1:DoFs
        for k=1:NumObe
            ADPR(k,i)= CoeModeS(k,i)^2/(NF(k,2)*2*pi); % (NF(k,2); 2nd NF
        end
    end
    ADPRDOF= sum(ADPR,1);
    ADPRDOF = normalize(ADPRDOF,'range');    

    CoeSecM =[];ADPRM=[]; 
    for i=1:NumSen
        CoeSecM = [CoeSecM CoeModeS(:,C(1,i))];  
        ADPRM = [ADPRM ADPRDOF(1,C(1,i))];
    end

    for i=1:NumSen
        CoeSelMMid = CoeSecM;
        CoeSelMMid(:,i) = [];
        ADPRMMid = ADPRM;
        ADPRMMid(:,i) = []; 
        
        [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
        MEdMid(1,i)= sum(CCC.^2)*sum(ADPRMMid); 
    end    
[~,IRep] = min(MEdMid);
MEdMidIni = MEdMid;
MEdMid(IRep) = [];   
SSCFSR = min(MEdMid);
Idel = find(MEdMidIni==SSCFSR,1, 'first');
end
