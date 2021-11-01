clear;
load Coe4ModesNorDam.mat X Y

MNum = 2; % Select the mode shape, 1:{1,1}; 2:{1,2}; 3{1,3}; 4{1,4};
CoeModeS = 10^6*X{1,MNum}; 

SNTol= 36;
CandiSL = 1:1:SNTol;

%% =========================== Load results ===============================

load OptResSSCGAFS.mat bestIndicesSSCFSGA optimalSSCFSGA
load OptResSSCGAFSR.mat bestIndicesSSCFSRGA optimalSSCFSRGA

%% *************************** Two Sensor locations ***********************
NumSen=2;

% GAR 
C =[];C = bestIndicesSSCFSGA{NumSen,1};
[optimalSSCFSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSRGA{NumSen,1};
[optimalSSCFSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y);

%% *************************** Three Sensor locations *********************
NumSen=3;

% GAR 
C =[];C = bestIndicesSSCFSGA{NumSen,1};
[optimalSSCFSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSRGA{NumSen,1};
[optimalSSCFSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y);

%% *************************** Four Sensor locations **********************
NumSen=4;

% GAR 
C =[];C = bestIndicesSSCFSGA{NumSen,1};
[optimalSSCFSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSRGA{NumSen,1};
[optimalSSCFSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y);

%% *************************** Five Sensor locations **********************
NumSen=5;

% GAR 
C =[];C = bestIndicesSSCFSGA{NumSen,1};
[optimalSSCFSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSRGA{NumSen,1};
[optimalSSCFSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y);

%% *************************** Six Sensor locations ***********************
NumSen=6;

% GAR 
C =[];C = bestIndicesSSCFSGA{NumSen,1};
[optimalSSCFSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSRGA{NumSen,1};
[optimalSSCFSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y);

%% *************************** Seven Sensor locations *********************
NumSen=7;

% GAR 
C =[];C = bestIndicesSSCFSGA{NumSen,1};
[optimalSSCFSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSRGA{NumSen,1};
[optimalSSCFSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y);

%% *************************** Eight Sensor locations *********************
NumSen=8;

% GAR 
C =[];C = bestIndicesSSCFSGA{NumSen,1};
[optimalSSCFSGA2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSRGA{NumSen,1};
[optimalSSCFSRGA2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y);


%% ======== Save the results  =============================================
bestIndicesSSCFSGA2=bestIndicesSSCFSGA;
for i = 2:8
    bestIndicesSSCFSGA2{i,1}(1,FSSenDel(1,i))= FSSenOptimal(1,i); 
    bestIndicesSSCFSGA2{i,1}=sort(bestIndicesSSCFSGA2{i,1},2);
end

bestIndicesSSCFSRGA2=bestIndicesSSCFSRGA;
for i = 2:8
    bestIndicesSSCFSRGA2{i,1}(1,FSRSenDel(1,i))= FSRSenOptimal(1,i); 
    bestIndicesSSCFSRGA2{i,1}=sort(bestIndicesSSCFSRGA2{i,1},2);
end

optimalSSCFSGA2 =optimalSSCFSGA2';
optimalSSCFSRGA2 =optimalSSCFSRGA2';

%%
save OptResSSCGAFS2.mat bestIndicesSSCFSGA2 optimalSSCFSGA2
save OptResSSCGAFSR2.mat bestIndicesSSCFSRGA2 optimalSSCFSRGA2

%% ======== Function 1  ===================================================
function [SSCMax,I,OptimalS,CandiSL0]= SSCMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y)
    CoeSelM =[]; 
    for i=1:NumSen
        CoeSelM = [CoeSelM CoeModeS(:,[C(1,i)])];  
    end

    for i=1:NumSen
        CoeSelMMid = CoeSelM;
        CoeSelMMid(:,i) = [];
        [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
        MEdMid(1,i)= sum(CCC.^2);   
    end    
    [SSCFS, I]= min(MEdMid);
    C(:,I)=[];
    CStable = C;

    DelN= size(CStable,2);
    for i =1:DelN
        CandiSL(CandiSL==C(1,i))=[];
    end
    CandiSL0 = CandiSL;
    
    NN = size(CandiSL,2);SSC =zeros(1,size(CandiSL,2));
    for i= 1:NN
        C = [CStable,CandiSL(1,i)];
        CoeSelM =[]; 
        for j=1:NumSen
            CoeSelM = [CoeSelM CoeModeS(:,[C(1,j)])];  
        end
        
        for j=1:NumSen
            CoeSelMMid = CoeSelM;
            CoeSelMMid(:,j) = [];
            [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
            MEdMid(1,j)= sum(CCC.^2); 
        end  
        
        if min(MEdMid)<SSCFS
           CandiSL0(:,CandiSL0==CandiSL(1,i))=[];
        else
           [~,~,CCC] = canoncorr(CoeSelM,Y);          
           SSC(1,i)=sum(CCC.^2);  
        end
    end 
    [SSCMax,Ind] = max(SSC);
    OptimalS = CandiSL(1,Ind);
end

%% ======== Function 2  ===================================================
function [SSCMax,Idel,OptimalS,CandiSL0]= SSCMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y)
    CoeSelM =[]; 
    for i=1:NumSen
        CoeSelM = [CoeSelM CoeModeS(:,[C(1,i)])];  
    end

    for i=1:NumSen
        CoeSelMMid = CoeSelM;
        CoeSelMMid(:,i) = [];
        [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
        MEdMid(1,i)= sum(CCC.^2);  
    end    
    [~,IRep] = min(MEdMid);
    MEdMidIni = MEdMid;
    MEdMid(:,IRep) = [];   
    SSCFSR = min(MEdMid);
    Idel= find(MEdMidIni==SSCFSR);
    C(:,Idel)=[];
    CStable = C;

    DelN= size(CStable,2);
    for i =1:DelN
        CandiSL(CandiSL==C(1,i))=[];
    end
    CandiSL0 = CandiSL;
    
    NN = size(CandiSL,2);SSC =zeros(1,size(CandiSL,2));
    for i= 1:NN
        C = [CStable,CandiSL(1,i)];
        CoeSelM =[]; 
        for j=1:NumSen
            CoeSelM = [CoeSelM CoeModeS(:,[C(1,j)])];  
        end
        
        for j=1:NumSen
            CoeSelMMid = CoeSelM;
            CoeSelMMid(:,j) = [];
            [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
            MEdMid(1,j)= sum(CCC.^2); 
        end  
        
        [~,IRep] = min(MEdMid);
        MEdMid(:,IRep) = [];   
        SSCFSR2 = min(MEdMid);
        if SSCFSR2<SSCFSR
           CandiSL0(:,CandiSL0==CandiSL(1,i))=[];
        else
           [~,~,CCC] = canoncorr(CoeSelM,Y);          
           SSC(1,i)=sum(CCC.^2);    
        end
    end 
    [SSCMax,Ind] = max(SSC);
    OptimalS = CandiSL(1,Ind);
end
