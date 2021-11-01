clear;
load Coe4Modes.mat Coe4Modes NF
Coe4ModeM =10^3*Coe4Modes{1,4};

SNTol= 36;
CandiSL = 1:1:SNTol;

%% =========================== Load results ===============================
load OptResFIMADPRGAFS.mat bestIndicesFIMADPR_FSGA optimalFIMADPR_FSGA
load OptResFIMADPRGAFSR.mat bestIndicesFIMADPR_FSRGA optimalFIMADPR_FSRGA

%% *************************** Four Sensor locations **********************
NumSen=4;

% GAR 
C =[];C = bestIndicesFIMADPR_FSGA{NumSen,1};
[optimalFIMADPR_FSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMADPRMax_FucFS(NumSen,C,Coe4ModeM,CandiSL,NF);

%% GAFSR
C =[];C = bestIndicesFIMADPR_FSRGA{NumSen,1};
[optimalFIMADPR_FSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMADPRMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL,NF);

%% *************************** Five Sensor locations **********************
NumSen=5;

% GAR 
C =[];C = bestIndicesFIMADPR_FSGA{NumSen,1};
[optimalFIMADPR_FSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMADPRMax_FucFS(NumSen,C,Coe4ModeM,CandiSL,NF);

%% GAFSR
C =[];C = bestIndicesFIMADPR_FSRGA{NumSen,1};
[optimalFIMADPR_FSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMADPRMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL,NF);

%% *************************** Six Sensor locations ***********************
NumSen=6;

% GAR 
C =[];C = bestIndicesFIMADPR_FSGA{NumSen,1};
[optimalFIMADPR_FSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMADPRMax_FucFS(NumSen,C,Coe4ModeM,CandiSL,NF);

%% GAFSR
C =[];C = bestIndicesFIMADPR_FSRGA{NumSen,1};
[optimalFIMADPR_FSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMADPRMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL,NF);

%% *************************** Seven Sensor locations *********************
NumSen=7;

% GAR 
C =[];C = bestIndicesFIMADPR_FSGA{NumSen,1};
[optimalFIMADPR_FSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMADPRMax_FucFS(NumSen,C,Coe4ModeM,CandiSL,NF);

%% GAFSR
C =[];C = bestIndicesFIMADPR_FSRGA{NumSen,1};
[optimalFIMADPR_FSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMADPRMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL,NF);

%% *************************** Eight Sensor locations *********************
NumSen=8;

% GAR 
C =[];C = bestIndicesFIMADPR_FSGA{NumSen,1};
[optimalFIMADPR_FSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMADPRMax_FucFS(NumSen,C,Coe4ModeM,CandiSL,NF);

%% GAFSR
C =[];C = bestIndicesFIMADPR_FSRGA{NumSen,1};
[optimalFIMADPR_FSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMADPRMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL,NF);

%% ======== Save the results  =============================================
bestIndicesFIMADPR_FSGA2=bestIndicesFIMADPR_FSGA;
for i = 4:8
    bestIndicesFIMADPR_FSGA2{i,1}(1,FSSenDel(1,i))= FSSenOptimal(1,i); 
    bestIndicesFIMADPR_FSGA2{i,1}=sort(bestIndicesFIMADPR_FSGA2{i,1},2);
end

bestIndicesFIMADPR_FSRGA2=bestIndicesFIMADPR_FSRGA;
for i = 4:8
    bestIndicesFIMADPR_FSRGA2{i,1}(1,FSRSenDel(1,i))= FSRSenOptimal(1,i); 
    bestIndicesFIMADPR_FSRGA2{i,1}=sort(bestIndicesFIMADPR_FSRGA2{i,1},2);
end

optimalFIMADPR_FSGA2 =optimalFIMADPR_FSGA2';
optimalFIMADPR_FSRGA2 =optimalFIMADPR_FSRGA2';

%%
save OptResFIMADPRGAFS2.mat bestIndicesFIMADPR_FSGA2 optimalFIMADPR_FSGA2
save OptResFIMADPRGAFSR2.mat bestIndicesFIMADPR_FSRGA2 optimalFIMADPR_FSRGA2

%% ======== Function 1  ===================================================
function [FIMADPRDetMax,I,OptimalS,CandiSL0]= FIMADPRMax_FucFS(NumSen,C,Coe4ModeM,CandiSL,NF)
    DOFs = size(Coe4ModeM,1);
    MSNum = size(NF,2);
    for i=1:DOFs
        for k=1:MSNum
            ADPR(i,k)= Coe4ModeM(i,k)^2/(NF(1,k)*2*pi);
        end
    end
    ADPRDOF= sum(ADPR,2);
    ADPRDOF = normalize(ADPRDOF,'range');    

    CoeSecM =[]; ADPRM = [];
    for i=1:NumSen
        CoeSecM = [CoeSecM;Coe4ModeM([C(1,i)],:)];  
        ADPRM = [ADPRM;ADPRDOF(C(1,i),1)];
    end

    for i=1:NumSen
        CoeSecMMid = CoeSecM;
        CoeSecMMid(i,:) = [];
        ADPRMMid = ADPRM;
        ADPRMMid(i,:) =[];
        MEdMid(1,i)= det(CoeSecMMid.'*CoeSecMMid)*sum(ADPRMMid); 
    end    
    [FIMADPRDetFS, I]= min(MEdMid);
    C(:,I)=[];
    CStable = C;

    DelN= size(CStable,2);
    for i =1:DelN
        CandiSL(CandiSL==C(1,i))=[];
    end
    CandiSL0 = CandiSL;
    
    NN = size(CandiSL,2);FIMADPRDet =zeros(1,size(CandiSL,2));
    for i= 1:NN
        C = [CStable,CandiSL(1,i)];
        CoeSecM =[]; ADPRM = []; 
        for j=1:NumSen
            CoeSecM = [CoeSecM;Coe4ModeM([C(1,j)],:)];  
            ADPRM = [ADPRM;ADPRDOF(C(1,j),1)];
        end
        
        for j=1:NumSen
            CoeSecMMid = CoeSecM;
            CoeSecMMid(j,:) = [];
            ADPRMMid = ADPRM;
            ADPRMMid(j,:) =[];
            MEdMid(1,j)= det(CoeSecMMid.'*CoeSecMMid)*sum(ADPRMMid); 
        end  
        
        if min(MEdMid)<FIMADPRDetFS
           CandiSL0(:,CandiSL0==CandiSL(1,i))=[];
        else
           FIMADPRDet(1,i)= det(CoeSecM.'*CoeSecM)*sum(ADPRM);    
        end
    end 
    [FIMADPRDetMax,Ind] = max(FIMADPRDet);
    OptimalS = CandiSL(1,Ind);
end

%% ======== Function 2  ===================================================
function [FIMADPRDetMax,Idel,OptimalS,CandiSL0]= FIMADPRMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL,NF)
    DOFs = size(Coe4ModeM,1);
    MSNum = size(NF,2);
    for i=1:DOFs
        for k=1:MSNum
            ADPR(i,k)= Coe4ModeM(i,k)^2/(NF(1,k)*2*pi);
        end
    end
    ADPRDOF= sum(ADPR,2);
    ADPRDOF = normalize(ADPRDOF,'range');    

    CoeSecM =[]; ADPRM = [];
    for i=1:NumSen
        CoeSecM = [CoeSecM;Coe4ModeM([C(1,i)],:)];  
        ADPRM = [ADPRM;ADPRDOF(C(1,i),1)];
    end

    for i=1:NumSen
        CoeSecMMid = CoeSecM;
        CoeSecMMid(i,:) = [];
        ADPRMMid = ADPRM;
        ADPRMMid(i,:) =[];
        MEdMid(1,i)= det(CoeSecMMid.'*CoeSecMMid)*sum(ADPRMMid); 
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
        CoeSecM =[]; ADPRM = []; 
        for j=1:NumSen
            CoeSecM = [CoeSecM;Coe4ModeM([C(1,j)],:)];  
            ADPRM = [ADPRM;ADPRDOF(C(1,j),1)];
        end
        
        for j=1:NumSen
            CoeSecMMid = CoeSecM;
            CoeSecMMid(j,:) = [];
            ADPRMMid = ADPRM;
            ADPRMMid(j,:) =[];
            MEdMid(1,j)= det(CoeSecMMid.'*CoeSecMMid)*sum(ADPRMMid); 
        end  
        
        [~,IRep] = min(MEdMid);
        MEdMid(:,IRep) = [];   
        FIMDetFSR2 = min(MEdMid);
        if FIMDetFSR2<FIMDetFSR
           CandiSL0(:,CandiSL0==CandiSL(1,i))=[];
        else
           FIMDet(1,i)= det(CoeSecM.'*CoeSecM)*sum(ADPRM);    
        end
    end 
    [FIMADPRDetMax,Ind] = max(FIMDet);
    OptimalS = CandiSL(1,Ind);
end
