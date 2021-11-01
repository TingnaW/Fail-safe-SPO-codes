clear;
load Coe4Modes.mat
Coe4ModeM = 10^3*Coe4Modes{1,4}; % Select the mode shapes corresponding to one temperature.
                            % 0:{1,1}; 5:{1,2}; 10{1,3}; 
                            % 15{1,4}; 20{1,5}; 25{1,6};
SN = size(Coe4ModeM,1);
CandidateSN = 1:1:SN;

%% ================ Select 4 sensors ======================================
NumSen=4;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[FIMDetFSEH{1,NumSen},LoopTFIMFSEH(NumSen,1)] = ExhaustiveS_FIMDetFSFuc(NumSen,C,CombN,Coe4ModeM);

[bestf, bestidx] = max(FIMDetFSEH{1,NumSen});    
optimalFIMFSEH(NumSen,1) = bestf
bestIndicesFIMFSEH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,FIMDetFSEH{1,NumSen})
xlabel('Combination number')
ylabel('Determinant of FIM')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(FIMDetFSEH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ================ Select 5 sensors ======================================
NumSen=5;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[FIMDetFSEH{1,NumSen},LoopTFIMFSEH(NumSen,1)] = ExhaustiveS_FIMDetFSFuc(NumSen,C,CombN,Coe4ModeM);

[bestf, bestidx] = max(FIMDetFSEH{1,NumSen});    
optimalFIMFSEH(NumSen,1) = bestf
bestIndicesFIMFSEH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,FIMDetFSEH{1,NumSen})
xlabel('Combination number')
ylabel('Determinant of FIM')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(FIMDetFSEH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ================ Select 6 sensors ======================================
NumSen=6;
C=[]; C = nchoosek(CandidateSN,NumSen);  CombN = size(C,1);

[FIMDetFSEH{1,NumSen},LoopTFIMFSEH(NumSen,1)] = ExhaustiveS_FIMDetFSFuc(NumSen,C,CombN,Coe4ModeM);

[bestf, bestidx] = max(FIMDetFSEH{1,NumSen});    
optimalFIMFSEH(NumSen,1) = bestf
bestIndicesFIMFSEH{NumSen,1} = C(bestidx,:)

AxisX=[];
AxisX = 1:CombN;
figure(NumSen+30)
plot(AxisX,FIMDetFSEH{1,NumSen})
xlabel('Combination number')
ylabel('Determinant of FIM')
set(gca,'FontSize',12)

[MaxK, Indall] = maxk(FIMDetFSEH{1,NumSen},10);
UniVal{1,NumSen} = unique(MaxK);
Ind = MaxK==UniVal{1,NumSen}(end,1);
SameLoc{1,NumSen}= C(Indall(Ind),:);

%% ======== Save the results  =============================================
save('OptResFIMExhFS','optimalFIMFSEH','bestIndicesFIMFSEH','LoopTFIMFSEH','FIMDetFSEH');

%% ======== Function  =====================================================
function [MEd,TimeP] = ExhaustiveS_FIMDetFSFuc(NumSen,C,CombN,Coe4ModeM)
tic
f = waitbar(0,sprintf('%7.4f%% iteration',0));
    for j = 1:CombN
        CoeSecM =[];
        for i=1:NumSen
            CoeSecM = [CoeSecM;Coe4ModeM([C(j,i)],:)];  
        end
        
        for i=1:NumSen
            CoeSecMMid = CoeSecM;
            CoeSecMMid(i,:) = [];
            MEdMid(1,i)= det(CoeSecMMid.'*CoeSecMMid);
        end  

        MEd(j,1)= min(MEdMid);            
        waitbar(j/CombN,f,sprintf('%7.4f%% iteration',j/CombN*100));
    end
close(f)
TimeP = toc;
end