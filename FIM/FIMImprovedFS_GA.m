clear;
load Coe4Modes.mat Coe4Modes NF
Coe4ModeM =10^3*Coe4Modes{1,4};

SNTol= 36;
CandiSL = 1:1:SNTol;

%% =========================== Load results ===============================

load OptResFIMGAFS.mat bestIndicesFIMFSGA optimalFIMFSGA
load OptResFIMGAFSR.mat bestIndicesFIMFSRGA optimalFIMFSRGA

%% *************************** Four Sensor locations **********************
NumSen=4;

% GAFS 
C =[];C = bestIndicesFIMFSGA{NumSen,1};
[optimalFIMFSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMMax_FucFS(NumSen,C,Coe4ModeM,CandiSL);

% GAFSR 
C =[];C = bestIndicesFIMFSRGA{NumSen,1};
[optimalFIMFSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL);

%% *************************** Five Sensor locations **********************
NumSen=5;

% GAFS 
C =[];C = bestIndicesFIMFSGA{NumSen,1};
[optimalFIMFSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMMax_FucFS(NumSen,C,Coe4ModeM,CandiSL);

% GAFSR 
C =[];C = bestIndicesFIMFSRGA{NumSen,1};
[optimalFIMFSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL);

%% *************************** Six Sensor locations ***********************
NumSen=6;

% GAFS 
C =[];C = bestIndicesFIMFSGA{NumSen,1};
[optimalFIMFSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMMax_FucFS(NumSen,C,Coe4ModeM,CandiSL);

% GAFSR 
C =[];C = bestIndicesFIMFSRGA{NumSen,1};
[optimalFIMFSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL);

%% *************************** Seven Sensor locations *********************
NumSen=7;

% GAFS 
C =[];C = bestIndicesFIMFSGA{NumSen,1};
[optimalFIMFSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMMax_FucFS(NumSen,C,Coe4ModeM,CandiSL);

% GAFSR 
C =[];C = bestIndicesFIMFSRGA{NumSen,1};
[optimalFIMFSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL);

%% *************************** Eight Sensor locations *********************
NumSen=8;

% GAFS 
C =[];C = bestIndicesFIMFSGA{NumSen,1};
[optimalFIMFSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMMax_FucFS(NumSen,C,Coe4ModeM,CandiSL);

% GAFSR 
C =[];C = bestIndicesFIMFSRGA{NumSen,1};
[optimalFIMFSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL);

%% ======== Save the results  =============================================
bestIndicesFIMFSGA2=bestIndicesFIMFSGA;
for i = 4:8
    bestIndicesFIMFSGA2{i,1}(1,FSSenDel(1,i))= FSSenOptimal(1,i); 
    bestIndicesFIMFSGA2{i,1}=sort(bestIndicesFIMFSGA2{i,1},2);
end

bestIndicesFIMFSRGA2=bestIndicesFIMFSRGA;
for i = 4:8
    bestIndicesFIMFSRGA2{i,1}(1,FSRSenDel(1,i))= FSRSenOptimal(1,i); 
    bestIndicesFIMFSRGA2{i,1}=sort(bestIndicesFIMFSRGA2{i,1},2);
end

optimalFIMFSGA2 =optimalFIMFSGA2';
optimalFIMFSRGA2 =optimalFIMFSRGA2';

%%
save OptResFIMGAFS2.mat bestIndicesFIMFSGA2 optimalFIMFSGA2
save OptResFIMGAFSR2.mat bestIndicesFIMFSRGA2 optimalFIMFSRGA2

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