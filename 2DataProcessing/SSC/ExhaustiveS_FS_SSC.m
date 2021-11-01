clear;
load Coe4ModesNorDam.mat X Y

MNum = 2; % Select the mode shape, 1:{1,1}; 2:{1,2}; 3{1,3}; 4{1,4};
CoeModeS = 10^6*X{1,MNum}; 
                                                                           
SN = size(CoeModeS,2);
CandidateSN = 1:1:SN;

%% ================ Select 2 sensors ======================================
NumSen=2;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCFSEH{1,NumSen},LoopTSSCFSEH(NumSen,1)] = ExhaustiveS_SSCFSFuc(NumSen,C,CombN,CoeModeS,Y);

[bestf, bestidx] = max(SSCFSEH{1,NumSen});    
optimalSSCFSEH(NumSen,1) = bestf
bestIndicesSSCFSEH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCFSEH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCFSEH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ================ Select 3 sensors ======================================
NumSen=3;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCFSEH{1,NumSen},LoopTSSCFSEH(NumSen,1)] = ExhaustiveS_SSCFSFuc(NumSen,C,CombN,CoeModeS,Y);

[bestf, bestidx] = max(SSCFSEH{1,NumSen});    
optimalSSCFSEH(NumSen,1) = bestf
bestIndicesSSCFSEH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCFSEH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCFSEH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ================ Select 4 sensors ======================================
NumSen=4;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCFSEH{1,NumSen},LoopTSSCFSEH(NumSen,1)] = ExhaustiveS_SSCFSFuc(NumSen,C,CombN,CoeModeS,Y);

[bestf, bestidx] = max(SSCFSEH{1,NumSen});    
optimalSSCFSEH(NumSen,1) = bestf
bestIndicesSSCFSEH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCFSEH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCFSEH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ================ Select 5 sensors ======================================
NumSen=5;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCFSEH{1,NumSen},LoopTSSCFSEH(NumSen,1)] = ExhaustiveS_SSCFSFuc(NumSen,C,CombN,CoeModeS,Y);

[bestf, bestidx] = max(SSCFSEH{1,NumSen});    
optimalSSCFSEH(NumSen,1) = bestf
bestIndicesSSCFSEH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCFSEH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCFSEH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ================ Select 6 sensors ======================================
NumSen=6;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCFSEH{1,NumSen},LoopTSSCFSEH(NumSen,1)] = ExhaustiveS_SSCFSFuc(NumSen,C,CombN,CoeModeS,Y);

[bestf, bestidx] = max(SSCFSEH{1,NumSen});    
optimalSSCFSEH(NumSen,1) = bestf
bestIndicesSSCFSEH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCFSEH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCFSEH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ======== Save the results  =============================================
save('OptResSSCExhFS','optimalSSCFSEH','bestIndicesSSCFSEH','LoopTSSCFSEH','SSCFSEH');

%% ======== Function  ===================================================== 
function [MEd,TimeP] = ExhaustiveS_SSCFSFuc(NumSen,C,CombN,CoeModeS,Y)
tic
f = waitbar(0,sprintf('%7.4f%% iteration',0));
    for j = 1:CombN
        CoeSelM =[];
        for i=1:NumSen
            CoeSelM = [CoeSelM CoeModeS(:,C(j,i))];  
        end
        
        for i=1:NumSen
            CoeSelMMid = CoeSelM;
            CoeSelMMid(:,i) = [];  
            
            [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
            MEdMid(1,i)= sum(CCC.^2);   
        end  

        MEd(j,1)= min(MEdMid);          
        waitbar(j/CombN,f,sprintf('%7.4f%% iteration',j/CombN*100));
    end
close(f)
TimeP = toc;
end