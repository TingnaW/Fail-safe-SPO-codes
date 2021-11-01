clear;
colors = lines(5);
FS = 15;
FS2 = 14;

%% =========================== Load results ===============================
% addpath('add the location of the following files')  
addpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/2DataProcessing/FIM')

load OptResFIMExh.mat bestIndicesFIMEH FIMDetEH optimalFIMEH
load OptResFIMExhFS2.mat bestIndicesFIMFSEH2 optimalFIMFSEH2
load OptResFIMExhFSR2.mat bestIndicesFIMFSREH2 optimalFIMFSREH2

load OptResFIMGA.mat bestIndicesFIMGA optimalFIMGA
load OptResFIMGAFS2.mat bestIndicesFIMFSGA2 optimalFIMFSGA2
load OptResFIMGAFSR2.mat bestIndicesFIMFSRGA2 optimalFIMFSRGA2

load Coe4Modes.mat Coe4Modes NF
Coe4ModeM =10^3*Coe4Modes{1,4};

% rmpath('add the location of these files')
rmpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/2DataProcessing/FIM')

%% *************************** Four Sensors *******************************
NumSen=4;
NumSen2 = NumSen-1;

figure(9)
histogram(FIMDetEH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Ehx
XLim =[optimalFIMEH(NumSen,1) optimalFIMEH(NumSen,1)];YLim = [0 4*10^4];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[optimalFIMFSEH2(NumSen,1) optimalFIMFSEH2(NumSen,1)];YLim = [0 4*10^4];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% GA
XLim =optimalFIMGA(NumSen,1);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFIMFSGA2(NumSen,1);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p1,p2,p4,p5],'Exhaustive search (ES)','Improved fail-safe ES','Genetic algorithm (GA)','Improved fail-safe GA','FontSize',FS2,'FontName','times')

%% *********************** After one sensor fail + Worse Case *************
% Exh
C =[];C = bestIndicesFIMEH{NumSen,1};
[FSFIMExh(1,NumSen), DelNumExh(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% ExhFS
C =[];C = bestIndicesFIMFSEH2{NumSen,1};
[FSFIMExhFS(1,NumSen), DelNumExhFS(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GA
C =[];C = bestIndicesFIMGA{NumSen,1};
[FSFIMGA(1,NumSen), DelNumGA(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFS 
C =[];C = bestIndicesFIMFSGA2{NumSen,1};
[FSFIMGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

figure(10) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
histogram(FIMDetEH{1,NumSen2},30,'FaceColor',colors(5,:)) 
hold on 

% Ehx
XLim =[FSFIMExh(1,NumSen) FSFIMExh(1,NumSen)];YLim = [0 6*10^3];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMExhFS(1,NumSen) FSFIMExhFS(1,NumSen)];YLim = [0 6*10^3];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% GA
XLim =FSFIMGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSFIMGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Five Sensors *******************************
NumSen=5;
NumSen2 = NumSen-1;

figure(1)
histogram(FIMDetEH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Ehx
XLim =[optimalFIMEH(NumSen,1) optimalFIMEH(NumSen,1)];YLim = [0 2*10^5];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[optimalFIMFSEH2(NumSen,1) optimalFIMFSEH2(NumSen,1)];YLim = [0 2*10^5];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% EhxFSR
XLim =[optimalFIMFSREH2(NumSen2,1) optimalFIMFSREH2(NumSen2,1)];YLim = [0 2*10^5];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalFIMGA(NumSen,1); YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFIMFSGA2(NumSen,1); YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =optimalFIMFSRGA2(NumSen2,1); YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p3,p6],'Improved fail-safe with redundancy ES','Improved fail-safe with redundancy GA','FontSize',FS2,'FontName','times')

%% *********************** After one sensor fail + Worse Case *************
% Exh
C =[];C = bestIndicesFIMEH{NumSen,1};
[FSFIMExh(1,NumSen), DelNumExh(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% ExhFS
C =[];C = bestIndicesFIMFSEH2{NumSen,1};
[FSFIMExhFS(1,NumSen), DelNumExhFS(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% ESFSR 
C =[];C = bestIndicesFIMFSREH2{NumSen2,1};
[FSFIMExhFSR(1,NumSen2), DelNumExhFSR(1,NumSen2),RepNumExhFSR(1,NumSen2)]= FIM_FucFSR(NumSen2,C,Coe4ModeM);

% GA
C =[];C = bestIndicesFIMGA{NumSen,1};
[FSFIMGA(1,NumSen), DelNumGA(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFS 
C =[];C = bestIndicesFIMFSGA2{NumSen,1};
[FSFIMGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFSR 
C =[];C = bestIndicesFIMFSRGA2{NumSen2,1};
[FSFIMGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= FIM_FucFSR(NumSen2,C,Coe4ModeM);

figure(2) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
histogram(FIMDetEH{1,NumSen2},30,'FaceColor',colors(5,:)) 
hold on 

% Ehx
XLim =[FSFIMExh(1,NumSen) FSFIMExh(1,NumSen)];YLim = [0 4*10^4];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMExhFS(1,NumSen) FSFIMExhFS(1,NumSen)];YLim = [0 4*10^4];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% EhxFSR
XLim =[FSFIMExhFSR(1,NumSen2) FSFIMExhFSR(1,NumSen2)];YLim = [0 4*10^4];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSFIMGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSFIMGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =FSFIMGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Six Sensors ********************************
NumSen=6;
NumSen2 = NumSen-1;

figure(3)
histogram(FIMDetEH{1,NumSen},30,'FaceColor',colors(5,:))
hold on 
% Ehx
XLim =[optimalFIMEH(NumSen,1) optimalFIMEH(NumSen,1)];YLim = [0 9*10^5];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[optimalFIMFSEH2(NumSen,1) optimalFIMFSEH2(NumSen,1)];YLim = [0 9*10^5];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% EhxFSR
XLim =[optimalFIMFSREH2(NumSen2,1) optimalFIMFSREH2(NumSen2,1)];YLim = [0 9*10^5];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalFIMGA(NumSen,1); YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFIMFSGA2(NumSen,1); YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =optimalFIMFSRGA2(NumSen2,1); YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')
ylim([0 9*10^5])

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case *************
% Exh
C =[];C = bestIndicesFIMEH{NumSen,1};
[FSFIMExh(1,NumSen), DelNumExh(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% ExhFS
C =[];C = bestIndicesFIMFSEH2{NumSen,1};
[FSFIMExhFS(1,NumSen), DelNumExhFS(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% ESFSR 
C =[];C = bestIndicesFIMFSREH2{NumSen2,1};
[FSFIMExhFSR(1,NumSen2), DelNumExhFSR(1,NumSen2),RepNumExhFSR(1,NumSen2)]= FIM_FucFSR(NumSen2,C,Coe4ModeM);

% GA
C =[];C = bestIndicesFIMGA{NumSen,1};
[FSFIMGA(1,NumSen), DelNumGA(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFS 
C =[];C = bestIndicesFIMFSGA2{NumSen,1};
[FSFIMGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFSR 
C =[];C = bestIndicesFIMFSRGA2{NumSen2,1};
[FSFIMGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= FIM_FucFSR(NumSen2,C,Coe4ModeM);

figure(4) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
histogram(FIMDetEH{1,NumSen2},30,'FaceColor',colors(5,:)) 
hold on 

% Ehx
XLim =[FSFIMExh(1,NumSen) FSFIMExh(1,NumSen)];YLim = [0 2*10^5];
p1=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMExhFS(1,NumSen) FSFIMExhFS(1,NumSen)];YLim = [0 2*10^5];
p2=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% EhxFSR
XLim =[FSFIMExhFSR(1,NumSen2) FSFIMExhFSR(1,NumSen2)];YLim = [0 2*10^5];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSFIMGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSFIMGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =FSFIMGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Seven Sensors ******************************
NumSen=7;
NumSen2 = NumSen-1;

figure(5)
hold on 

% EhxFSR
XLim =[optimalFIMFSREH2(NumSen2,1) optimalFIMFSREH2(NumSen2,1)];YLim = [0 2*10^5];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalFIMGA(NumSen,1); YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFIMFSGA2(NumSen,1); YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =optimalFIMFSRGA2(NumSen2,1); YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case *************
% ExhFSR 
C =[];C = bestIndicesFIMFSREH2{NumSen2,1};
[FSFIMExhFSR(1,NumSen2), DelNumExhFSR(1,NumSen2),RepNumExhFSR(1,NumSen2)]= FIM_FucFSR(NumSen2,C,Coe4ModeM);

% GA
C =[];C = bestIndicesFIMGA{NumSen,1};
[FSFIMGA(1,NumSen), DelNumGA(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFS 
C =[];C = bestIndicesFIMFSGA2{NumSen,1};
[FSFIMGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFSR 
C =[];C = bestIndicesFIMFSRGA2{NumSen2,1};
[FSFIMGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= FIM_FucFSR(NumSen2,C,Coe4ModeM);

figure(6) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
hold on 

% EhxFSR
XLim =[FSFIMExhFSR(1,NumSen2) FSFIMExhFSR(1,NumSen2)];YLim = [0 2*10^5];
p3=plot(XLim,YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSFIMGA(1,NumSen);YLim = 0;
p4=plot(XLim,YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =FSFIMGAFS(1,NumSen);YLim = 0;
p5=plot(XLim,YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =FSFIMGAFSR(1,NumSen2);YLim = 0;
p6=plot(XLim,YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('DFIM');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Eight Sensor locations *********************
NumSen=8;
NumSen2 = NumSen-1;

figure(7)
hold on 
% GA
XLim =[optimalFIMGA(NumSen,1) optimalFIMGA(NumSen,1)];YLim = [0 2*10^5];
p2=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% GAFS
XLim =[optimalFIMFSGA2(NumSen,1) optimalFIMFSGA2(NumSen,1)];YLim = [0 2*10^5];
p3=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% GAFSR
XLim =[optimalFIMFSRGA2(NumSen2,1) optimalFIMFSRGA2(NumSen2,1)];YLim = [0 2*10^5];
p4=plot(XLim,YLim,'--','Color',colors(3,:),'LineWidth',2);
xlabel('DFIM');ylabel('Frequency');set(gca,'FontSize',FS)

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p2,p3,p4],'Genetic algorithm','Improved Fail-safe (GA)','Improved Fail-safe with redundancy (GA)','FontSize',FS2)

%% *********************** After one sensor fail + Worse Case *************

% GA
C =[];C = bestIndicesFIMGA{NumSen,1};
[FSFIMGA(1,NumSen), DelNumGA(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFS 
C =[];C = bestIndicesFIMFSGA2{NumSen,1};
[FSFIMGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFSR 
C =[];C = bestIndicesFIMFSRGA2{NumSen2,1};
[FSFIMGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= FIM_FucFSR(NumSen2,C,Coe4ModeM);

figure(8) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
hold on 

% GA
XLim =[FSFIMGA(1,NumSen) FSFIMGA(1,NumSen)];YLim = [0 2*10^5];
p2=plot(XLim,YLim,'-.','Color',colors(1,:),'LineWidth',2);

% GAFS
XLim =[FSFIMGAFS(1,NumSen) FSFIMGAFS(1,NumSen)];YLim = [0 2*10^5];
p3=plot(XLim,YLim,'--','Color',colors(2,:),'LineWidth',2);

% GAFSR
XLim =[FSFIMGAFSR(1,NumSen2) FSFIMGAFSR(1,NumSen2)];YLim = [0 2*10^5];
p4=plot(XLim,YLim,'--','Color',colors(3,:),'LineWidth',2);
xlabel('DFIM');ylabel('Frequency');set(gca,'FontSize',FS)

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p2,p3,p4],'Genetic algorithm','Improved Fail-safe (GA)','Improved Fail-safe with redundancy (GA)','FontSize',FS2)

%% *************************** Nine Sensor locations **********************
% GAFSR 
NumSen=9;
NumSen2 = NumSen-1;
C =[];C = bestIndicesFIMFSRGA2{NumSen2,1};
[FSFIMGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= FIM_FucFSR(NumSen2,C,Coe4ModeM);

%% ======== Save the results  =============================================
save FIMGASensorDistribution.mat bestIndicesFIMGA DelNumGA bestIndicesFIMFSGA2 DelNumGAFS bestIndicesFIMFSRGA2 DelNumGAFSR RepNumGAFSR

%% ======== Function 1  ===================================================
function [FIMDetFS, I]= FIM_FucFS(NumSen,C,Coe4ModeM)
CoeSecM =[]; 
    for i=1:NumSen
        CoeSecM = [CoeSecM;Coe4ModeM([C(1,i)],:)];  
    end

    for i=1:NumSen
        CoeSecMMid = CoeSecM;
        CoeSecMMid(i,:) = [];
        MEdMid(1,i)= det(CoeSecMMid.'*CoeSecMMid);
    end    
[FIMDetFS, I]= min(MEdMid);
end

%% ======== Function 2  ===================================================
function [FIMDetFSR, Idel, IRep]= FIM_FucFSR(NumSen,C,Coe4ModeM)
CoeSecM =[]; 
    for i=1:NumSen
        CoeSecM = [CoeSecM;Coe4ModeM([C(1,i)],:)];  
    end

    for i=1:NumSen
        CoeSecMMid = CoeSecM;
        CoeSecMMid(i,:) = [];
        MEdMid(1,i)= det(CoeSecMMid.'*CoeSecMMid);
    end    
[~,IRep] = min(MEdMid);
MEdMidIni = MEdMid;
MEdMid(IRep) = [];   
FIMDetFSR = min(MEdMid);
Idel = find(MEdMidIni==FIMDetFSR,1, 'first');
end

