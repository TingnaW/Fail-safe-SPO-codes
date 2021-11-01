clear;
load Coe4Modes.mat Coe4Modes NF
Coe4ModeM =10^3*Coe4Modes{1,4};

SNTol= 36;
CandiSL = 1:1:SNTol;

%% =========================== Load results ===============================

load OptResFIMExhFS bestIndicesFIMFSEH optimalFIMFSEH  
load OptResFIMExhFSR bestIndicesFIMFSREH optimalFIMFSREH 

%% *************************** Four Sensor locations **********************
NumSen=4;

% GAR 
C =[];C = bestIndicesFIMFSEH{NumSen,1};
[optimalFIMFSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMMax_FucFS(NumSen,C,Coe4ModeM,CandiSL);

% GAFSR 
C =[];C = bestIndicesFIMFSREH{NumSen,1};
[optimalFIMFSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL);

%% *************************** Five Sensor locations **********************
NumSen=5;

% GAR 
C =[];C = bestIndicesFIMFSEH{NumSen,1};
[optimalFIMFSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMMax_FucFS(NumSen,C,Coe4ModeM,CandiSL);

% GAFSR 
C =[];C = bestIndicesFIMFSREH{NumSen,1};
[optimalFIMFSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL);

%% *************************** Six Sensor locations ***********************
NumSen=6;

% GAR 
C =[];C = bestIndicesFIMFSEH{NumSen,1};
[optimalFIMFSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMMax_FucFS(NumSen,C,Coe4ModeM,CandiSL);

% GAFSR 
C =[];C = bestIndicesFIMFSREH{NumSen,1};
[optimalFIMFSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL);

%% ======== Save the results  =============================================
bestIndicesFIMFSEH2=bestIndicesFIMFSEH;
for i = 4:6
    bestIndicesFIMFSEH2{i,1}(1,FSSenDel(1,i))= FSSenOptimal(1,i); 
    bestIndicesFIMFSEH2{i,1}=sort(bestIndicesFIMFSEH2{i,1},2);
end

bestIndicesFIMFSREH2=bestIndicesFIMFSREH;
for i = 4:6
    bestIndicesFIMFSREH2{i,1}(1,FSRSenDel(1,i))= FSRSenOptimal(1,i); 
    bestIndicesFIMFSREH2{i,1}=sort(bestIndicesFIMFSREH2{i,1},2);
end

optimalFIMFSEH2 =optimalFIMFSEH2';
optimalFIMFSREH2 =optimalFIMFSREH2';

%%
save OptResFIMExhFS2.mat bestIndicesFIMFSEH2 optimalFIMFSEH2
save OptResFIMExhFSR2.mat bestIndicesFIMFSREH2 optimalFIMFSREH2

%% ======== Function 1  ===================================================
function [FIMDetMax,I,OptimalS,CandiSL0]= FIMMax_FucFS(NumSen,C,Coe4ModeM,CandiSL)
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
    C(:,I)=[];
    CStable = C;

    DelN= size(CStable,2);
    for i =1:DelN
        CandiSL(CandiSL==C(1,i))=[];
    end
    CandiSL0 = CandiSL;
    
    NN = size(CandiSL,2);FIMDet =zeros(1,size(CandiSL,2));
    for i= 1:NN
        C = [CStable,CandiSL(1,i)];
        CoeSecM =[]; 
        for j=1:NumSen
            CoeSecM = [CoeSecM;Coe4ModeM([C(1,j)],:)];  
        end
        
        for j=1:NumSen
            CoeSecMMid = CoeSecM;
            CoeSecMMid(j,:) = [];
            MEdMid(1,j)= det(CoeSecMMid.'*CoeSecMMid);
        end  
        
        if min(MEdMid)<FIMDetFS
           CandiSL0(:,CandiSL0==CandiSL(1,i))=[];
        else
           FIMDet(1,i)= det(CoeSecM.'*CoeSecM);    
        end
    end 
    [FIMDetMax,Ind] = max(FIMDet);
    OptimalS = CandiSL(1,Ind);
end

%% ======== Function 2  ===================================================
function [FIMDetMax,Idel,OptimalS,CandiSL0]= FIMMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL)
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
    MEdMid(:,IRep) = [];   
    FIMDetFSR = min(MEdMid);
    Idel= find(MEdMidIni==FIMDetFSR);
    C(:,Idel)=[];
    CStable = C;

    DelN= size(CStable,2);
    for i =1:DelN
        CandiSL(CandiSL==C(1,i))=[];
    end
    CandiSL0 = CandiSL;
    
    NN = size(CandiSL,2);FIMDet =zeros(1,size(CandiSL,2));
    for i= 1:NN
        C = [CStable,CandiSL(1,i)];
        CoeSecM =[]; 
        for j=1:NumSen
            CoeSecM = [CoeSecM;Coe4ModeM([C(1,j)],:)];  
        end
        
        for j=1:NumSen
            CoeSecMMid = CoeSecM;
            CoeSecMMid(j,:) = [];
            MEdMid(1,j)= det(CoeSecMMid.'*CoeSecMMid);
        end  
        
        [~,IRep] = min(MEdMid);
        MEdMid(:,IRep) = [];   
        FIMDetFSR2 = min(MEdMid);
        if FIMDetFSR2<FIMDetFSR
           CandiSL0(:,CandiSL0==CandiSL(1,i))=[];
        else
           FIMDet(1,i)= det(CoeSecM.'*CoeSecM);    
        end
    end 
    [FIMDetMax,Ind] = max(FIMDet);
    OptimalS = CandiSL(1,Ind);
end
