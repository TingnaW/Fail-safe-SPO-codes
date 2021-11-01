clear;
load Coe4ModesNorDam.mat X Y

MNum = 2; % Select the mode shape, 1:{1,1}; 2:{1,2}; 3{1,3}; 4{1,4};
CoeModeS = 10^6*X{1,MNum}; 
                                                                           
SN = size(CoeModeS,2);
CandidateSN = 1:1:SN;

%% ================ Select 2 sensor locations =============================
NumSen=2;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCFSREH{1,NumSen},LoopTSSCFSREH(NumSen,1)] = ExhaustiveS_SSCFSRFuc(NumSen,C,CombN,CoeModeS,Y);

[bestf, bestidx] = max(SSCFSREH{1,NumSen});    
optimalSSCFSREH(NumSen,1) = bestf
bestIndicesSSCFSREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCFSREH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCFSREH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ================ Select 3 sensor locations =============================
NumSen=3;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCFSREH{1,NumSen},LoopTSSCFSREH(NumSen,1)] = ExhaustiveS_SSCFSRFuc(NumSen,C,CombN,CoeModeS,Y);

[bestf, bestidx] = max(SSCFSREH{1,NumSen});    
optimalSSCFSREH(NumSen,1) = bestf
bestIndicesSSCFSREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCFSREH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCFSREH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ================ Select 4 sensor locations =============================
NumSen=4;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCFSREH{1,NumSen},LoopTSSCFSREH(NumSen,1)] = ExhaustiveS_SSCFSRFuc(NumSen,C,CombN,CoeModeS,Y);

[bestf, bestidx] = max(SSCFSREH{1,NumSen});    
optimalSSCFSREH(NumSen,1) = bestf
bestIndicesSSCFSREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCFSREH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCFSREH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ================ Select 5 sensor locations =============================
NumSen=5;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCFSREH{1,NumSen},LoopTSSCFSREH(NumSen,1)] = ExhaustiveS_SSCFSRFuc(NumSen,C,CombN,CoeModeS,Y);

[bestf, bestidx] = max(SSCFSREH{1,NumSen});    
optimalSSCFSREH(NumSen,1) = bestf
bestIndicesSSCFSREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCFSREH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCFSREH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ================ Select 6 sensor locations =============================
NumSen=6;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[SSCFSREH{1,NumSen},LoopTSSCFSREH(NumSen,1)] = ExhaustiveS_SSCFSRFuc(NumSen,C,CombN,CoeModeS,Y);

[bestf, bestidx] = max(SSCFSREH{1,NumSen});    
optimalSSCFSREH(NumSen,1) = bestf
bestIndicesSSCFSREH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,SSCFSREH{1,NumSen})
xlabel('Combination number')
ylabel('SSC')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(SSCFSREH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ======== Save the results  =============================================
save('OptResSSCExhFSR','optimalSSCFSREH','bestIndicesSSCFSREH','LoopTSSCFSREH','SSCFSREH');

%% ======== Function  ===================================================== 
function [MEd,TimeP] = ExhaustiveS_SSCFSRFuc(NumSen,C,CombN,CoeModeS,Y)
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
        [~,Ind] = min(MEdMid);
        MEdMid(Ind) = [];
        MEd(j,1)= min(MEdMid);          
        waitbar(j/CombN,f,sprintf('%7.4f%% iteration',j/CombN*100));
    end
close(f)
TimeP = toc;
end