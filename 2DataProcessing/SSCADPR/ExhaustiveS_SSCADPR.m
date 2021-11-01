clear;
load Coe4ModesNorDam.mat

MNum = 2; % Select the mode shape, 1:{1,1}; 2:{1,2}; 3{1,3}; 4{1,4};
CoeModeS = 10^6*X{1,MNum}; 
                                                                           
SN = size(CoeModeS,2);
CandidateSN = 1:1:SN;

%% ================ Select 1 sensors ======================================
NumSen=1;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCADPREH{1,NumSen},LoopTSSCADPREH(NumSen,1)] = ExhaustiveS_SSCADPRFuc(NumSen,C,CombN,CoeModeS,Y,NF);

[bestf, bestidx] = max(SSCADPREH{1,NumSen});    
optimalSSCADPREH(NumSen,1) = bestf
bestIndicesSSCADPREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+20)
plot(AxisX,SSCADPREH{1,NumSen})
xlabel('Combination number')
ylabel('Effective independence index')
set(gca,'FontSize',12)

%% ================ Select 2 sensors ======================================
NumSen=2;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCADPREH{1,NumSen},LoopTSSCADPREH(NumSen,1)] = ExhaustiveS_SSCADPRFuc(NumSen,C,CombN,CoeModeS,Y,NF);

[bestf, bestidx] = max(SSCADPREH{1,NumSen});    
optimalSSCADPREH(NumSen,1) = bestf
bestIndicesSSCADPREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+20)
plot(AxisX,SSCADPREH{1,NumSen})
xlabel('Combination number')
ylabel('Effective independence index')
set(gca,'FontSize',12)

%% ================ Select 3 sensors ======================================
NumSen=3;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCADPREH{1,NumSen},LoopTSSCADPREH(NumSen,1)] = ExhaustiveS_SSCADPRFuc(NumSen,C,CombN,CoeModeS,Y,NF);

[bestf, bestidx] = max(SSCADPREH{1,NumSen});    
optimalSSCADPREH(NumSen,1) = bestf
bestIndicesSSCADPREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+20)
plot(AxisX,SSCADPREH{1,NumSen})
xlabel('Combination number')
ylabel('Effective independence index')
set(gca,'FontSize',12)

%% ================ Select 4 sensors ======================================
NumSen=4;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCADPREH{1,NumSen},LoopTSSCADPREH(NumSen,1)] = ExhaustiveS_SSCADPRFuc(NumSen,C,CombN,CoeModeS,Y,NF);

[bestf, bestidx] = max(SSCADPREH{1,NumSen});    
optimalSSCADPREH(NumSen,1) = bestf
bestIndicesSSCADPREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+20)
plot(AxisX,SSCADPREH{1,NumSen})
xlabel('Combination number')
ylabel('Effective independence index')
set(gca,'FontSize',12)

%% ================ Select 5 sensors ======================================
NumSen=5;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCADPREH{1,NumSen},LoopTSSCADPREH(NumSen,1)] = ExhaustiveS_SSCADPRFuc(NumSen,C,CombN,CoeModeS,Y,NF);

[bestf, bestidx] = max(SSCADPREH{1,NumSen});    
optimalSSCADPREH(NumSen,1) = bestf
bestIndicesSSCADPREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+20)
plot(AxisX,SSCADPREH{1,NumSen})
xlabel('Combination number')
ylabel('Effective independence index')
set(gca,'FontSize',12)

%% ================ Select 6 sensors ======================================
NumSen=6;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCADPREH{1,NumSen},LoopTSSCADPREH(NumSen,1)] = ExhaustiveS_SSCADPRFuc(NumSen,C,CombN,CoeModeS,Y,NF);

[bestf, bestidx] = max(SSCADPREH{1,NumSen});    
optimalSSCADPREH(NumSen,1) = bestf
bestIndicesSSCADPREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+20)
plot(AxisX,SSCADPREH{1,NumSen})
xlabel('Combination number')
ylabel('Effective independence index')
set(gca,'FontSize',12)

%% ======== Save the results  =============================================
save('OptResSSCADPRExh','optimalSSCADPREH','bestIndicesSSCADPREH','LoopTSSCADPREH','SSCADPREH');

%% ======== Function  ===================================================== 
function [CCCSqSum,TimeP] = ExhaustiveS_SSCADPRFuc(NumSen,C,CombN,CoeModeS,Y,NF)
tic
f = waitbar(0,sprintf('%7.4f%% iteration',0));
    NumObe = size(CoeModeS,1);
    DoFs = size(CoeModeS,2);
    for i=1:DoFs
        for k=1:NumObe
            ADPR(k,i)= CoeModeS(k,i)^2/(NF(k,2)*2*pi); % (NF(k,2); 2nd NF
        end
    end
    ADPRDOF= sum(ADPR,1);
    ADPRDOF = normalize(ADPRDOF,'range');        


    for j = 1:CombN
        CoeSelM =[];ADPRM=[];
        for i=1:NumSen
            CoeSelM = [CoeSelM CoeModeS(:,C(j,i))];  
            ADPRM = [ADPRM ADPRDOF(1,C(j,i))];
        end

        [~,~,CCC] = canoncorr(CoeSelM,Y);           
        CCCSqSum(j,1)=sum(CCC.^2)*sum(ADPRM);
        waitbar(j/CombN,f,sprintf('%7.4f%% iteration',j/CombN*100));
    end
close(f)
TimeP = toc;
end