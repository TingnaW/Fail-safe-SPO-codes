clear;
load Coe4Modes.mat Coe4Modes NF
Coe4ModeM =10^3*Coe4Modes{1,4};

SNTol= 36;
CandiSL = 1:1:SNTol;

%% =========================== Load results ===============================
load OptResFIMADPRExhFS bestIndicesFIMADPR_FSEH optimalFIMADPR_FSEH  
load OptResFIMADPRExhFSR bestIndicesFIMADPR_FSREH optimalFIMADPR_FSREH 

%% *************************** Four Sensor locations **********************
NumSen=4;

% GAR 
C =[];C = bestIndicesFIMADPR_FSEH{NumSen,1};
[optimalFIMADPR_FSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMADPRMax_FucFS(NumSen,C,Coe4ModeM,CandiSL,NF);

%% GAFSR
C =[];C = bestIndicesFIMADPR_FSREH{NumSen,1};
[optimalFIMADPR_FSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMADPRMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL,NF);

%% *************************** Five Sensor locations **********************
NumSen=5;

% GAR 
C =[];C = bestIndicesFIMADPR_FSEH{NumSen,1};
[optimalFIMADPR_FSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMADPRMax_FucFS(NumSen,C,Coe4ModeM,CandiSL,NF);

%% GAFSR
C =[];C = bestIndicesFIMADPR_FSREH{NumSen,1};
[optimalFIMADPR_FSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMADPRMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL,NF);

%% *************************** Six Sensor locations ***********************
NumSen=6;

% GAR 
C =[];C = bestIndicesFIMADPR_FSEH{NumSen,1};
[optimalFIMADPR_FSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= FIMADPRMax_FucFS(NumSen,C,Coe4ModeM,CandiSL,NF);

%% GAFSR
C =[];C = bestIndicesFIMADPR_FSREH{NumSen,1};
[optimalFIMADPR_FSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= FIMADPRMax_FucFSR(NumSen,C,Coe4ModeM,CandiSL,NF);

%% ======== Save the results  =============================================
bestIndicesFIMADPR_FSEH2=bestIndicesFIMADPR_FSEH;
for i = 4:6
    bestIndicesFIMADPR_FSEH2{i,1}(1,FSSenDel(1,i))= FSSenOptimal(1,i); 
    bestIndicesFIMADPR_FSEH2{i,1}=sort(bestIndicesFIMADPR_FSEH2{i,1},2);
end

bestIndicesFIMADPR_FSREH2=bestIndicesFIMADPR_FSREH;
for i = 4:6
    bestIndicesFIMADPR_FSREH2{i,1}(1,FSRSenDel(1,i))= FSRSenOptimal(1,i); 
    bestIndicesFIMADPR_FSREH2{i,1}=sort(bestIndicesFIMADPR_FSREH2{i,1},2);
end

optimalFIMADPR_FSEH2 =optimalFIMADPR_FSEH2';
optimalFIMADPR_FSREH2 =optimalFIMADPR_FSREH2';

%%
save OptResFIMADPRExhFS2.mat bestIndicesFIMADPR_FSEH2 optimalFIMADPR_FSEH2
save OptResFIMADPRExhFSR2.mat bestIndicesFIMADPR_FSREH2 optimalFIMADPR_FSREH2

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
