clear;
load Coe4Modes.mat Coe4Modes NF
Coe4ModeM =10^3*Coe4Modes{1,4};

colors = lines(5);
FS = 15;
FS2 = 14;

%% =========================== Part 2 FIM - ADPR =================================
load OptResFIMADPRExh.mat bestIndicesFIMADPR FIMADPR optimalFIMADPR
load OptResFIMADPRExhFS2.mat bestIndicesFIMADPRFSEH2 optimalFSFIMADPREH2
load OptResFIMADPRExhFSR2.mat bestIndicesFIMADPRFSREH2 optimalFSRFIMADPREH2

load OptResFIMADPRGA.mat bestIndicesFIMADPRGA optimalFIMADPRGA
load OptResFIMADPRGAFS2.mat bestIndicesFIMADPRFS2 optimalFSFIMADPR
load OptResFIMADPRGAFSR2.mat bestIndicesFIMADPRFSR2 optimalFSRFIMADPR

%% *************************** Four Sensors *******************************
NumSen=4;
NumSen2 = NumSen-1;

figure(9)
histogram(log(FIMADPR{1,NumSen}),30,'FaceColor',colors(5,:))
hold on 
% Ehx
XLim =[optimalFIMADPR(NumSen,1) optimalFIMADPR(NumSen,1)];YLim = [0 7*10^3];
p1=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[optimalFSFIMADPREH2(NumSen,1) optimalFSFIMADPREH2(NumSen,1)];YLim = [0 7*10^3];
p2=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);

% GA
XLim =optimalFIMADPRGA(NumSen,1);YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFSFIMADPR(NumSen,1);YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM weighted by ADPR)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
% legend([p1,p2,p4,p5],'Exhaustive search (ES)','Improved fail-safe ES','Genetic algorithm (GA)','Improved fail-safe GA','FontSize',FS2,'FontName','times')

%% *********************** After one sensor fail + Worse Case ***************
% Exh
C =[];C = bestIndicesFIMADPR{NumSen,1};
[FSFIMADPRExh(1,NumSen), DelNumADPRExh(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% Exh
C =[];C = bestIndicesFIMADPRFSEH2{NumSen,1};
[FSFIMADPRExhFS(1,NumSen), DelNumADPRExhFS(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% GA
C =[];C = bestIndicesFIMADPRGA{NumSen,1};
[FSFIMADPRGA(1,NumSen), DelNumADPRGA(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% GAFS 
C =[];C = bestIndicesFIMADPRFS2{NumSen,1};
[FSFIMADPRGAFS(1,NumSen), DelNumADPRGAFS(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

figure(10) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
histogram(log(FIMADPR{1,NumSen2}),30,'FaceColor',colors(5,:)) % ^^^^^^^^^^^^^^^^
hold on 

% Ehx
XLim =[FSFIMADPRExh(1,NumSen) FSFIMADPRExh(1,NumSen)];YLim = [0 10^3];
p1=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMADPRExhFS(1,NumSen) FSFIMADPRExhFS(1,NumSen)];YLim = [0 10^3];
p2=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);

% GA
XLim =FSFIMADPRGA(1,NumSen);YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSFIMADPRGAFS(1,NumSen);YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM weighted by ADPR)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Five Sensors *******************************
NumSen=5;
NumSen2 = NumSen-1;

figure(1)
histogram(log(FIMADPR{1,NumSen}),30,'FaceColor',colors(5,:))
hold on 
% Ehx
XLim =[optimalFIMADPR(NumSen,1) optimalFIMADPR(NumSen,1)];YLim = [0 5*10^4];
p1=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[optimalFSFIMADPREH2(NumSen,1) optimalFSFIMADPREH2(NumSen,1)];YLim = [0 5*10^4];
p2=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[optimalFSRFIMADPREH2(NumSen2,1) optimalFSRFIMADPREH2(NumSen2,1)];YLim = [0 5*10^4];
p3=plot(log(XLim),YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalFIMADPRGA(NumSen,1);YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFSFIMADPR(NumSen,1);YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =optimalFSRFIMADPR(NumSen2,1);YLim = 0;
p6=plot(log(XLim),YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM weighted by ADPR)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
% legend([p3,p6],'Improved fail-safe with redundancy ES','Improved fail-safe with redundancy GA','FontSize',FS2,'FontName','times')

%% *********************** After one sensor fail + Worse Case ***************
% Exh
C =[];C = bestIndicesFIMADPR{NumSen,1};
[FSFIMADPRExh(1,NumSen), DelNumADPRExh(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% Exh
C =[];C = bestIndicesFIMADPRFSEH2{NumSen,1};
[FSFIMADPRExhFS(1,NumSen), DelNumADPRExhFS(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% ExhFSR
C =[];C = bestIndicesFIMADPRFSREH2{NumSen2,1};
[FSFIMADPRExhFSR(1,NumSen2), DelNumADPRExhFSR(1,NumSen2),RepNumADPRExhFSR(1,NumSen2)]= FIMADPRFucFSR(NumSen2,C,Coe4ModeM,NF);

% GA
C =[];C = bestIndicesFIMADPRGA{NumSen,1};
[FSFIMADPRGA(1,NumSen), DelNumADPRGA(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% GAFS 
C =[];C = bestIndicesFIMADPRFS2{NumSen,1};
[FSFIMADPRGAFS(1,NumSen), DelNumADPRGAFS(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% GAFSR
C =[];C = bestIndicesFIMADPRFSR2{NumSen2,1};
[FSFIMADPRGAFSR(1,NumSen2), DelNumADPRGAFSR(1,NumSen2),RepNumADPRGAFSR(1,NumSen2)]= FIMADPRFucFSR(NumSen2,C,Coe4ModeM,NF);

figure(2) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
histogram(log(FIMADPR{1,NumSen2}),30,'FaceColor',colors(5,:)) % ^^^^^^^^^^^^^^^^
hold on 

% Ehx
XLim =[FSFIMADPRExh(1,NumSen) FSFIMADPRExh(1,NumSen)];YLim = [0 7*10^3];
p1=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMADPRExhFS(1,NumSen) FSFIMADPRExhFS(1,NumSen)];YLim = [0 7*10^3];
p2=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMADPRExhFSR(1,NumSen2) FSFIMADPRExhFSR(1,NumSen2)];YLim = [0 7*10^3];
p3=plot(log(XLim),YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSFIMADPRGA(1,NumSen);YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSFIMADPRGAFS(1,NumSen);YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =FSFIMADPRGAFSR(1,NumSen2);YLim = 0;
p6=plot(log(XLim),YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM weighted by ADPR)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Six Sensors *******************************
NumSen=6;
NumSen2 = NumSen-1;

figure(3)
histogram(log(FIMADPR{1,NumSen}),30,'FaceColor',colors(5,:))
hold on 
% Ehx
XLim =[optimalFIMADPR(NumSen,1) optimalFIMADPR(NumSen,1)];YLim = [0 2.7*10^5];
p1=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[optimalFSFIMADPREH2(NumSen,1) optimalFSFIMADPREH2(NumSen,1)];YLim = [0 2.7*10^5];
p2=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);

% ExhFSR
XLim =[optimalFSRFIMADPREH2(NumSen2,1) optimalFSRFIMADPREH2(NumSen2,1)];YLim = [0 2.7*10^5];
p3=plot(log(XLim),YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalFIMADPRGA(NumSen,1);YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFSFIMADPR(NumSen,1);YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =optimalFSRFIMADPR(NumSen2,1);YLim = 0;
p6=plot(log(XLim),YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM weighted by ADPR)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')
ylim([0 2.7*10^5])

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case ***************
% Exh
C =[];C = bestIndicesFIMADPR{NumSen,1};
[FSFIMADPRExh(1,NumSen), DelNumADPRExh(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% Exh
C =[];C = bestIndicesFIMADPRFSEH2{NumSen,1};
[FSFIMADPRExhFS(1,NumSen), DelNumADPRExhFS(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% ExhFSR
C =[];C = bestIndicesFIMADPRFSREH2{NumSen2,1};
[FSFIMADPRExhFSR(1,NumSen2), DelNumADPRExhFSR(1,NumSen2),RepNumADPRExhFSR(1,NumSen2)]= FIMADPRFucFSR(NumSen2,C,Coe4ModeM,NF);

% GA
C =[];C = bestIndicesFIMADPRGA{NumSen,1};
[FSFIMADPRGA(1,NumSen), DelNumADPRGA(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% GAFS 
C =[];C = bestIndicesFIMADPRFS2{NumSen,1};
[FSFIMADPRGAFS(1,NumSen), DelNumADPRGAFS(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% GAFSR
C =[];C = bestIndicesFIMADPRFSR2{NumSen2,1};
[FSFIMADPRGAFSR(1,NumSen2), DelNumADPRGAFSR(1,NumSen2),RepNumADPRGAFSR(1,NumSen2)]= FIMADPRFucFSR(NumSen2,C,Coe4ModeM,NF);

figure(4) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
histogram(log(FIMADPR{1,NumSen2}),30,'FaceColor',colors(5,:)) % ^^^^^^^^^^^^^^^^
hold on 

% Ehx
XLim =[FSFIMADPRExh(1,NumSen) FSFIMADPRExh(1,NumSen)];YLim = [0 5*10^4];
p1=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMADPRExhFS(1,NumSen) FSFIMADPRExhFS(1,NumSen)];YLim = [0 5*10^4];
p2=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMADPRExhFSR(1,NumSen2) FSFIMADPRExhFSR(1,NumSen2)];YLim = [0 5*10^4];
p3=plot(log(XLim),YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSFIMADPRGA(1,NumSen);YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSFIMADPRGAFS(1,NumSen);YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =FSFIMADPRGAFSR(1,NumSen2);YLim = 0;
p6=plot(log(XLim),YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM weighted by ADPR)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Seven Sensors *******************************
NumSen=7;
NumSen2 = NumSen-1;

figure(5)
hold on 

% ExhFSR
XLim =[optimalFSRFIMADPREH2(NumSen2,1) optimalFSRFIMADPREH2(NumSen2,1)];YLim = [0 3*10^5];
p3=plot(log(XLim),YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalFIMADPRGA(NumSen,1);YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFSFIMADPR(NumSen,1);YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =optimalFSRFIMADPR(NumSen2,1);YLim = 0;
p6=plot(log(XLim),YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM weighted by ADPR)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case ***************
% ExhFSR
C =[];C = bestIndicesFIMADPRFSREH2{NumSen2,1};
[FSFIMADPRExhFSR(1,NumSen2), DelNumADPRExhFSR(1,NumSen2),RepNumADPRExhFSR(1,NumSen2)]= FIMADPRFucFSR(NumSen2,C,Coe4ModeM,NF);

% GA
C =[];C = bestIndicesFIMADPRGA{NumSen,1};
[FSFIMADPRGA(1,NumSen), DelNumADPRGA(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% GAFS 
C =[];C = bestIndicesFIMADPRFS2{NumSen,1};
[FSFIMADPRGAFS(1,NumSen), DelNumADPRGAFS(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% GAFSR
C =[];C = bestIndicesFIMADPRFSR2{NumSen2,1};
[FSFIMADPRGAFSR(1,NumSen2), DelNumADPRGAFSR(1,NumSen2),RepNumADPRGAFSR(1,NumSen2)]= FIMADPRFucFSR(NumSen2,C,Coe4ModeM,NF);

figure(6) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
hold on 

% EhxFS
XLim =[FSFIMADPRExhFSR(1,NumSen2) FSFIMADPRExhFSR(1,NumSen2)];YLim = [0 3*10^5];
p3=plot(log(XLim),YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSFIMADPRGA(1,NumSen);YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSFIMADPRGAFS(1,NumSen);YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =FSFIMADPRGAFSR(1,NumSen2);YLim = 0;
p6=plot(log(XLim),YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM weighted by ADPR)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Eight Sensors *******************************
NumSen=8;
NumSen2 = NumSen-1;

figure(7)
hold on 

% GA
XLim =[optimalFIMADPRGA(NumSen,1) optimalFIMADPRGA(NumSen,1)];YLim = [0 3*10^5];
p2=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% GAR: Fail-safe
XLim =[optimalFSFIMADPR(NumSen,1) optimalFSFIMADPR(NumSen,1)];YLim = [0 3*10^5];
p3=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);
xlabel('log(DFIM weighted by ADPR)');ylabel('Frequency');set(gca,'FontSize',FS)

% GAR: Fail-safe with redundancy
XLim =[optimalFSRFIMADPR(NumSen2,1) optimalFSRFIMADPR(NumSen2,1)];YLim = [0 3*10^5];
p4=plot(log(XLim),YLim,'--','Color',colors(3,:),'LineWidth',2);
xlabel('log(DFIM weighted by ADPR)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p2,p3,p4],'Genetic algorithm','Improved Fail-safe (GA)','Improved Fail-safe with redundancy (GA)','FontSize',FS2)

%% *********************** After one sensor fail + Worse Case ***************

% GA
C =[];C = bestIndicesFIMADPRGA{NumSen,1};
[FSFIMADPRGA(1,NumSen), DelNumADPRGA(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% GAR 
C =[];C = bestIndicesFIMADPRFS2{NumSen,1};
[FSFIMADPRGAFS(1,NumSen), DelNumADPRGAFS(1,NumSen)]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF);

% GAFSR
C =[];C = bestIndicesFIMADPRFSR2{NumSen2,1};
[FSFIMADPRGAFSR(1,NumSen2), DelNumADPRGAFSR(1,NumSen2),RepNumADPRGAFSR(1,NumSen2)]= FIMADPRFucFSR(NumSen2,C,Coe4ModeM,NF);

figure(8) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
hold on 

% GA
XLim =[FSFIMADPRGA(1,NumSen) FSFIMADPRGA(1,NumSen)];YLim = [0 3*10^5];
p2=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);
% GAR: Fail-safe
XLim =[FSFIMADPRGAFS(1,NumSen) FSFIMADPRGAFS(1,NumSen)];YLim = [0 3*10^5];
p3=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);
xlabel('log(DFIM weighted by ADPR)');ylabel('Frequency');set(gca,'FontSize',FS)

% GAR: Fail-safe with redundancy
XLim =[FSFIMADPRGAFSR(1,NumSen2) FSFIMADPRGAFSR(1,NumSen2)];YLim = [0 3*10^5];
p4=plot(log(XLim),YLim,'--','Color',colors(3,:),'LineWidth',2);
xlabel('log(DFIM weighted by ADPR)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p2,p3,p4],'Genetic algorithm','Improved Fail-safe (GA)','Improved Fail-safe with redundancy (GA)','FontSize',FS2)

%%
% GAFSR
NumSen=9;
NumSen2 = NumSen-1;
C =[];C = bestIndicesFIMADPRFSR2{NumSen2,1};
[FSFIMADPRGAFSR(1,NumSen2), DelNumADPRGAFSR(1,NumSen2),RepNumADPRGAFSR(1,NumSen2)]= FIMADPRFucFSR(NumSen2,C,Coe4ModeM,NF);

%% 
% save FIMADPRGASensorDistribution.mat bestIndicesFIMADPRGA DelNumADPRGA bestIndicesFIMADPRFS2 DelNumADPRGAFS bestIndicesFIMADPRFSR2 DelNumADPRGAFSR RepNumADPRGAFSR

%%
Ans1= optimalFSFIMADPREH2-FSFIMADPRExhFS';
Ans2= optimalFSRFIMADPREH2 -FSFIMADPRExhFSR';
close all
figure(1)
plot(Ans1)
hold on
plot(Ans2)
hold off

Ans3= optimalFSFIMADPR-FSFIMADPRGAFS';
Ans4= optimalFSRFIMADPR -FSFIMADPRGAFSR';
figure(2)
plot(Ans3)
hold on
plot(Ans4)
hold off

%%
function [FIMADPRR, I]= FIMADPRFucFS(NumSen,C,Coe4ModeM,NF)

    DOFs = size(Coe4ModeM,1);
    MSNum = size(NF,2);
    ADPR=[];               % **********************************************
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

%%
function [FIMADPRR, Idel, IRep]= FIMADPRFucFSR(NumSen,C,Coe4ModeM,NF)

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
