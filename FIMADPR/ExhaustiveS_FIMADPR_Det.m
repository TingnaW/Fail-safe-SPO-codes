clear;
load Coe4Modes.mat
Coe4ModeM = 10^3*Coe4Modes{1,4}; % Select the mode shapes corresponding to one temperature.
                            % 0:{1,1}; 5:{1,2}; 10{1,3}; 
                            % 15{1,4}; 20{1,5}; 25{1,6};
SN = size(Coe4ModeM,1);
CandidateSN = 1:1:SN;

%% ================ Select 3 sensors ======================================
NumSen=3;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[FIMADPREH{1,NumSen},LoopTFIMADPREH(NumSen,1)] = ExhaustiveS_FIMADPRFuc(NumSen,C,CombN,Coe4ModeM,NF);

[bestf, bestidx] = max(FIMADPREH{1,NumSen});    
optimalFIMADPREH(NumSen,1) = bestf
bestIndicesFIMADPREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,FIMADPREH{1,NumSen})
xlabel('Combination number')
ylabel('Effective independence index')
set(gca,'FontSize',12)

%% ================ Select 4 sensors ======================================
NumSen=4;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[FIMADPREH{1,NumSen},LoopTFIMADPREH(NumSen,1)] = ExhaustiveS_FIMADPRFuc(NumSen,C,CombN,Coe4ModeM,NF);

[bestf, bestidx] = max(FIMADPREH{1,NumSen});    
optimalFIMADPREH(NumSen,1) = bestf
bestIndicesFIMADPREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,FIMADPREH{1,NumSen})
xlabel('Combination number')
ylabel('Effective independence index')
set(gca,'FontSize',12)

%% ================ Select 5 sensors ======================================
NumSen=5;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[FIMADPREH{1,NumSen},LoopTFIMADPREH(NumSen,1)] = ExhaustiveS_FIMADPRFuc(NumSen,C,CombN,Coe4ModeM,NF);

[bestf, bestidx] = max(FIMADPREH{1,NumSen});    
optimalFIMADPREH(NumSen,1) = bestf
bestIndicesFIMADPREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,FIMADPREH{1,NumSen})
xlabel('Combination number')
ylabel('Effective independence index')
set(gca,'FontSize',12)

%% ================ Select 6 sensors ======================================
NumSen=6;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[FIMADPREH{1,NumSen},LoopTFIMADPREH(NumSen,1)] = ExhaustiveS_FIMADPRFuc(NumSen,C,CombN,Coe4ModeM,NF);

[bestf, bestidx] = max(FIMADPREH{1,NumSen});    
optimalFIMADPREH(NumSen,1) = bestf
bestIndicesFIMADPREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,FIMADPREH{1,NumSen})
xlabel('Combination number')
ylabel('Effective independence index')
set(gca,'FontSize',12)

%% ======== Save the results  =============================================
save('OptResFIMADPRExh','optimalFIMADPREH','bestIndicesFIMADPREH','LoopTFIMADPREH','FIMADPREH');

%% ======== Function  =====================================================
function [MEd,TimeP] = ExhaustiveS_FIMADPRFuc(NumSen,C,CombN,Coe4ModeM,NF)
tic
f = waitbar(0,sprintf('%7.4f%% iteration',0));
    DOFs = size(Coe4ModeM,1);
    MSNum = size(Coe4ModeM,2);
    for i=1:DOFs
        for k=1:MSNum
            ADPR(i,k)= Coe4ModeM(i,k)^2/(NF(1,k)*2*pi);
        end
    end
    ADPRDOF= sum(ADPR,2);
    ADPRDOF = normalize(ADPRDOF,'range');        

    for j = 1:CombN
        CoeSecM =[]; ADPRM =[]; 
        for i=1:NumSen
            CoeSecM = [CoeSecM;Coe4ModeM(C(j,i),:)]; 
            ADPRM = [ADPRM;ADPRDOF(C(j,i),1)];
        end
        
        MEd(j,1)= det(CoeSecM.'*CoeSecM)*sum(ADPRM);                          
        waitbar(j/CombN,f,sprintf('%7.4f%% iteration',j/CombN*100));
    end
close(f)
TimeP = toc;
end