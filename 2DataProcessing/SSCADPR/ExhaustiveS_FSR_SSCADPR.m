clear;
load Coe4ModesNorDam.mat

MNum = 2; % Select the mode shape, 1:{1,1}; 2:{1,2}; 3{1,3}; 4{1,4};
CoeModeS = 10^6*X{1,MNum}; 
                                                                           
SN = size(CoeModeS,2);
CandidateSN = 1:1:SN;

%% ================ Select 2 sensors ======================================
NumSen=2;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCADPRFSREH{1,NumSen},LoopTSSCADPRFSREH(NumSen,1)] = ExhaustiveS_SSCADPRFSRFuc(NumSen,C,CombN,CoeModeS,Y,NF);

[bestf, bestidx] = max(SSCADPRFSREH{1,NumSen});    
optimalSSCADPRFSREH(NumSen,1) = bestf
bestIndicesSSCADPRFSREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCADPRFSREH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCADPRFSREH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ================ Select 3 sensors ======================================
NumSen=3;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCADPRFSREH{1,NumSen},LoopTSSCADPRFSREH(NumSen,1)] = ExhaustiveS_SSCADPRFSRFuc(NumSen,C,CombN,CoeModeS,Y,NF);

[bestf, bestidx] = max(SSCADPRFSREH{1,NumSen});    
optimalSSCADPRFSREH(NumSen,1) = bestf
bestIndicesSSCADPRFSREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCADPRFSREH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCADPRFSREH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ================ Select 4 sensors ======================================
NumSen=4;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCADPRFSREH{1,NumSen},LoopTSSCADPRFSREH(NumSen,1)] = ExhaustiveS_SSCADPRFSRFuc(NumSen,C,CombN,CoeModeS,Y,NF);

[bestf, bestidx] = max(SSCADPRFSREH{1,NumSen});    
optimalSSCADPRFSREH(NumSen,1) = bestf
bestIndicesSSCADPRFSREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCADPRFSREH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCADPRFSREH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ================ Select 5 sensors ======================================
NumSen=5;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCADPRFSREH{1,NumSen},LoopTSSCADPRFSREH(NumSen,1)] = ExhaustiveS_SSCADPRFSRFuc(NumSen,C,CombN,CoeModeS,Y,NF);

[bestf, bestidx] = max(SSCADPRFSREH{1,NumSen});    
optimalSSCADPRFSREH(NumSen,1) = bestf
bestIndicesSSCADPRFSREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCADPRFSREH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCADPRFSREH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ================ Select 6 sensors ======================================
NumSen=6;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCADPRFSREH{1,NumSen},LoopTSSCADPRFSREH(NumSen,1)] = ExhaustiveS_SSCADPRFSRFuc(NumSen,C,CombN,CoeModeS,Y,NF);

[bestf, bestidx] = max(SSCADPRFSREH{1,NumSen});    
optimalSSCADPRFSREH(NumSen,1) = bestf
bestIndicesSSCADPRFSREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCADPRFSREH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCADPRFSREH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ======== Save the results  =============================================
save('OptResSSCADPRExhFSR','optimalSSCADPRFSREH','bestIndicesSSCADPRFSREH','LoopTSSCADPRFSREH','SSCADPRFSREH');

%% ======== Function  ===================================================== 
function [MEd,TimeP] = ExhaustiveS_SSCADPRFSRFuc(NumSen,C,CombN,CoeModeS,Y,NF)
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
        
        for i=1:NumSen
            CoeSelMMid = CoeSelM;
            CoeSelMMid(:,i) = [];  
            ADPRMMid = ADPRM;                    
            ADPRMMid(:,i) = [];                            
            
            [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
            MEdMid(1,i)= sum(CCC.^2)*sum(ADPRMMid);    
        end  
        [~,Ind] = min(MEdMid);
        MEdMid(Ind) = [];
        MEd(j,1)= min(MEdMid);          
        waitbar(j/CombN,f,sprintf('%7.4f%% iteration',j/CombN*100));
    end
close(f)
TimeP = toc;
end