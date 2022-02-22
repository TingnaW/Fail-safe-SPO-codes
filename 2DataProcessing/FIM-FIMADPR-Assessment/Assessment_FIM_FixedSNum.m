clear;
load Coe4Modes.mat Coe4Modes NF
Coe4ModeM =10^3*Coe4Modes{1,4};

colors = lines(5);
FS = 15;
FS2 = 14;
%% =========================== Part 1 FIM =================================
load OptResFIMExh.mat bestIndicesFIM FIMDet optimalFIMEH
load OptResFIMExhFS2.mat bestIndicesFIMFSEH2 optimalFSFIMEH2
load OptResFIMExhFSR2.mat bestIndicesFIMFSREH2 optimalFSRFIMEH2

load OptResFIMGA.mat bestIndicesFIMGA optimalFIMGA
load OptResFIMGAFS2.mat bestIndicesFIMFS2 optimalFSFIM
load OptResFIMGAFSR2.mat bestIndicesFIMFSR2 optimalFSRFIM

%% *************************** Four Sensors *******************************
NumSen=4;
NumSen2 = NumSen-1;

figure(9)
histogram(log(FIMDet{1,NumSen}),30,'FaceColor',colors(5,:))
hold on 
% Ehx
XLim =[optimalFIMEH(NumSen,1) optimalFIMEH(NumSen,1)];YLim = [0 10^4];
p1=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[optimalFSFIMEH2(NumSen,1) optimalFSFIMEH2(NumSen,1)];YLim = [0 10^4];
p2=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);

% GA
XLim =optimalFIMGA(NumSen,1);YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe
XLim =optimalFSFIM(NumSen,1);YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p1,p2,p4,p5],'Exhaustive search (ES)','Improved fail-safe ES','Genetic algorithm (GA)','Improved fail-safe GA','FontSize',FS2,'FontName','times', 'Location','northwest')

%% *********************** After one sensor fail + Worse Case ***************
% Exh
C =[];C = bestIndicesFIM{NumSen,1};
[FSFIMExh(1,NumSen), DelNumExh(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% ExhFS
C =[];C = bestIndicesFIMFSEH2{NumSen,1};
[FSFIMExhFS(1,NumSen), DelNumExhFS(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GA
C =[];C = bestIndicesFIMGA{NumSen,1};
[FSFIMGA(1,NumSen), DelNumGA(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFS 
C =[];C = bestIndicesFIMFS2{NumSen,1};
[FSFIMGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

figure(10) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
histogram(log(FIMDet{1,NumSen2}),30,'FaceColor',colors(5,:)) % ^^^^^^^^^^^^^^^^
hold on 

% Ehx
XLim =[FSFIMExh(1,NumSen) FSFIMExh(1,NumSen)];YLim = [0 1.3*10^3];
p1=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMExhFS(1,NumSen) FSFIMExhFS(1,NumSen)];YLim = [0 1.3*10^3];
p2=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);

% GA
XLim =FSFIMGA(1,NumSen);YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe
XLim =FSFIMGAFS(1,NumSen);YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')
ylim([0 1.3*10^3])

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Five Sensors *******************************
NumSen=5;
NumSen2 = NumSen-1;

figure(1)
histogram(log(FIMDet{1,NumSen}),30,'FaceColor',colors(5,:))
hold on 
% Ehx
XLim =[optimalFIMEH(NumSen,1) optimalFIMEH(NumSen,1)];YLim = [0 6*10^4];
p1=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[optimalFSFIMEH2(NumSen,1) optimalFSFIMEH2(NumSen,1)];YLim = [0 6*10^4];
p2=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);

% EhxFSR
XLim =[optimalFSRFIMEH2(NumSen2,1) optimalFSRFIMEH2(NumSen2,1)];YLim = [0 6*10^4];
p3=plot(log(XLim),YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalFIMGA(NumSen,1); YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFSFIM(NumSen,1); YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =optimalFSRFIM(NumSen2,1); YLim = 0;
p6=plot(log(XLim),YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p3,p6],'Improved fail-safe with redundancy ES','Improved fail-safe with redundancy GA','FontSize',FS2,'FontName','times', 'Location','northwest')

%% *********************** After one sensor fail + Worse Case ***************
% Exh
C =[];C = bestIndicesFIM{NumSen,1};
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
C =[];C = bestIndicesFIMFS2{NumSen,1};
[FSFIMGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFSR 
C =[];C = bestIndicesFIMFSR2{NumSen2,1};
[FSFIMGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= FIM_FucFSR(NumSen2,C,Coe4ModeM);

figure(2) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
histogram(log(FIMDet{1,NumSen2}),30,'FaceColor',colors(5,:)) % ^^^^^^^^^^^^^^^^
hold on 

% Ehx
XLim =[FSFIMExh(1,NumSen) FSFIMExh(1,NumSen)];YLim = [0 0.9*10^4];
p1=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMExhFS(1,NumSen) FSFIMExhFS(1,NumSen)];YLim = [0 0.9*10^4];
p2=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);

% EhxFSR
XLim =[FSFIMExhFSR(1,NumSen2) FSFIMExhFSR(1,NumSen2)];YLim = [0 0.9*10^4];
p3=plot(log(XLim),YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSFIMGA(1,NumSen);YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe
XLim =FSFIMGAFS(1,NumSen);YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =FSFIMGAFSR(1,NumSen2);YLim = 0;
p6=plot(log(XLim),YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')
ylim([0 0.9*10^4])

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Six Sensors *******************************
NumSen=6;
NumSen2 = NumSen-1;

figure(3)
histogram(log(FIMDet{1,NumSen}),30,'FaceColor',colors(5,:))
hold on 
% Ehx
XLim =[optimalFIMEH(NumSen,1) optimalFIMEH(NumSen,1)];YLim = [0 2.7*10^5];
p1=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[optimalFSFIMEH2(NumSen,1) optimalFSFIMEH2(NumSen,1)];YLim = [0 2.7*10^5];
p2=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);

% EhxFSR
XLim =[optimalFSRFIMEH2(NumSen2,1) optimalFSRFIMEH2(NumSen2,1)];YLim = [0 2.7*10^5];
p3=plot(log(XLim),YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalFIMGA(NumSen,1); YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFSFIM(NumSen,1); YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =optimalFSRFIM(NumSen2,1); YLim = 0;
p6=plot(log(XLim),YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')
ylim([0 2.7*10^5])

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case ***************
% Exh
C =[];C = bestIndicesFIM{NumSen,1};
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
C =[];C = bestIndicesFIMFS2{NumSen,1};
[FSFIMGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFSR 
C =[];C = bestIndicesFIMFSR2{NumSen2,1};
[FSFIMGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= FIM_FucFSR(NumSen2,C,Coe4ModeM);

figure(4) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
histogram(log(FIMDet{1,NumSen2}),30,'FaceColor',colors(5,:)) % ^^^^^^^^^^^^^^^^
hold on 

% Ehx
XLim =[FSFIMExh(1,NumSen) FSFIMExh(1,NumSen)];YLim = [0 5*10^4];
p1=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% EhxFS
XLim =[FSFIMExhFS(1,NumSen) FSFIMExhFS(1,NumSen)];YLim = [0 5*10^4];
p2=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);

% EhxFSR
XLim =[FSFIMExhFSR(1,NumSen2) FSFIMExhFSR(1,NumSen2)];YLim = [0 5*10^4];
p3=plot(log(XLim),YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSFIMGA(1,NumSen);YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe
XLim =FSFIMGAFS(1,NumSen);YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =FSFIMGAFSR(1,NumSen2);YLim = 0;
p6=plot(log(XLim),YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Seven Sensors *******************************
NumSen=7;
NumSen2 = NumSen-1;

figure(5)
hold on 

% EhxFSR
XLim =[optimalFSRFIMEH2(NumSen2,1) optimalFSRFIMEH2(NumSen2,1)];YLim = [0 2*10^5];
p3=plot(log(XLim),YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =optimalFIMGA(NumSen,1); YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAFS
XLim =optimalFSFIM(NumSen,1); YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =optimalFSRFIM(NumSen2,1); YLim = 0;
p6=plot(log(XLim),YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *********************** After one sensor fail + Worse Case ***************
% ESFSR 
C =[];C = bestIndicesFIMFSREH2{NumSen2,1};
[FSFIMExhFSR(1,NumSen2), DelNumExhFSR(1,NumSen2),RepNumExhFSR(1,NumSen2)]= FIM_FucFSR(NumSen2,C,Coe4ModeM);

% GA
C =[];C = bestIndicesFIMGA{NumSen,1};
[FSFIMGA(1,NumSen), DelNumGA(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFS 
C =[];C = bestIndicesFIMFS2{NumSen,1};
[FSFIMGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFSR 
C =[];C = bestIndicesFIMFSR2{NumSen2,1};
[FSFIMGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= FIM_FucFSR(NumSen2,C,Coe4ModeM);

figure(6) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
hold on 

% EhxFSR
XLim =[FSFIMExhFSR(1,NumSen2) FSFIMExhFSR(1,NumSen2)];YLim = [0 2*10^5];
p3=plot(log(XLim),YLim,':','Color',colors(3,:),'LineWidth',2);

% GA
XLim =FSFIMGA(1,NumSen);YLim = 0;
p4=plot(log(XLim),YLim,'*','Color',colors(1,:),'MarkerSize',12,'LineWidth',2);

% GAR: Fail-safe
XLim =FSFIMGAFS(1,NumSen);YLim = 0;
p5=plot(log(XLim),YLim,'*','Color',colors(2,:),'MarkerSize',12,'LineWidth',2);

% GAFSR
XLim =FSFIMGAFSR(1,NumSen2);YLim = 0;
p6=plot(log(XLim),YLim,'*','Color',colors(3,:),'MarkerSize',12,'LineWidth',2);
xlabel('log(DFIM)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;

%% *************************** Eight Sensor locations *******************************
NumSen=8;
NumSen2 = NumSen-1;

figure(7)
hold on 
% GA
XLim =[optimalFIMGA(NumSen,1) optimalFIMGA(NumSen,1)];YLim = [0 2*10^5];
p2=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% GAR: Fail-safe
XLim =[optimalFSFIM(NumSen,1) optimalFSFIM(NumSen,1)];YLim = [0 2*10^5];
p3=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =[optimalFSRFIM(NumSen2,1) optimalFSRFIM(NumSen2,1)];YLim = [0 2*10^5];
p4=plot(log(XLim),YLim,'--','Color',colors(3,:),'LineWidth',2);
xlabel('log(DFIM)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p2,p3,p4],'Genetic algorithm','Improved Fail-safe (GA)','Improved Fail-safe with redundancy (GA)','FontSize',FS2)

%% *********************** After one sensor fail + Worse Case ***************

% GA
C =[];C = bestIndicesFIMGA{NumSen,1};
[FSFIMGA(1,NumSen), DelNumGA(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAR 
C =[];C = bestIndicesFIMFS2{NumSen,1};
[FSFIMGAFS(1,NumSen), DelNumGAFS(1,NumSen)]= FIM_FucFS(NumSen,C,Coe4ModeM);

% GAFSR 
C =[];C = bestIndicesFIMFSR2{NumSen2,1};
[FSFIMGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= FIM_FucFSR(NumSen2,C,Coe4ModeM);

figure(8) % ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
hold on 

% GA
XLim =[FSFIMGA(1,NumSen) FSFIMGA(1,NumSen)];YLim = [0 2*10^5];
p2=plot(log(XLim),YLim,'-.','Color',colors(1,:),'LineWidth',2);

% GAR: Fail-safe
XLim =[FSFIMGAFS(1,NumSen) FSFIMGAFS(1,NumSen)];YLim = [0 2*10^5];
p3=plot(log(XLim),YLim,'--','Color',colors(2,:),'LineWidth',2);

% GAR: Fail-safe with redundancy
XLim =[FSFIMGAFSR(1,NumSen2) FSFIMGAFSR(1,NumSen2)];YLim = [0 2*10^5];
p4=plot(log(XLim),YLim,'--','Color',colors(3,:),'LineWidth',2);
xlabel('log(DFIM)');ylabel('Frequency');set(gca,'FontSize',FS,'FontName','times')

ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p2,p3,p4],'Genetic algorithm','Improved Fail-safe (GA)','Improved Fail-safe with redundancy (GA)','FontSize',FS2)

%% 
% GAFSR 
NumSen=9;
NumSen2 = NumSen-1;
C =[];C = bestIndicesFIMFSR2{NumSen2,1};
[FSFIMGAFSR(1,NumSen2), DelNumGAFSR(1,NumSen2),RepNumGAFSR(1,NumSen2)]= FIM_FucFSR(NumSen2,C,Coe4ModeM);

%% 
% save FIMGASensorDistribution.mat bestIndicesFIMGA DelNumGA bestIndicesFIMFS2 DelNumGAFS bestIndicesFIMFSR2 DelNumGAFSR RepNumGAFSR

%%
Ans1= optimalFSFIMEH2-FSFIMExhFS';
Ans2= optimalFSRFIMEH2-FSFIMExhFSR';
close all
figure(1)
plot(Ans1)
hold on
plot(Ans2)
hold off

Ans3= optimalFSFIM-FSFIMGAFS';
Ans4= optimalFSRFIM-FSFIMGAFSR';
figure(2)
plot(Ans3)
hold on
plot(Ans4)
hold off

%%
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

%%
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

