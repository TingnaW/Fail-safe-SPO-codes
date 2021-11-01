clear;
load Coe4ModesNorDam.mat X Y

MNum = 2; % Select the mode shape, 1:{1,1}; 2:{1,2}; 3{1,3}; 4{1,4};
CoeModeS = 10^6*X{1,MNum}; 

SNTol= 36;
CandiSL = 1:1:SNTol;

%% =========================== Load results ===============================

load OptResSSCExhFS.mat bestIndicesSSCFSEH optimalSSCFSEH 
load OptResSSCExhFSR.mat bestIndicesSSCFSREH optimalSSCFSREH

%% *************************** Two Sensor locations ***********************
NumSen=2;

% GAR 
C =[];C = bestIndicesSSCFSEH{NumSen,1};
[optimalSSCFSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSREH{NumSen,1};
[optimalSSCFSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y);

%% *************************** Three Sensor locations *********************
NumSen=3;

% GAR 
C =[];C = bestIndicesSSCFSEH{NumSen,1};
[optimalSSCFSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSREH{NumSen,1};
[optimalSSCFSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y);

%% *************************** Four Sensor locations **********************
NumSen=4;

% GAR 
C =[];C = bestIndicesSSCFSEH{NumSen,1};
[optimalSSCFSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSREH{NumSen,1};
[optimalSSCFSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y);

%% *************************** Five Sensor locations **********************
NumSen=5;

% GAR 
C =[];C = bestIndicesSSCFSEH{NumSen,1};
[optimalSSCFSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSREH{NumSen,1};
[optimalSSCFSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y);

%% *************************** Six Sensor locations ***********************
NumSen=6;

% GAR 
C =[];C = bestIndicesSSCFSEH{NumSen,1};
[optimalSSCFSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y);

% GAFSR 
C =[];C = bestIndicesSSCFSREH{NumSen,1};
[optimalSSCFSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y);

%% ======== Save the results  =============================================
bestIndicesSSCFSEH2=bestIndicesSSCFSEH;
for i = 2:6
    bestIndicesSSCFSEH2{i,1}(1,FSSenDel(1,i))= FSSenOptimal(1,i); 
    bestIndicesSSCFSEH2{i,1}=sort(bestIndicesSSCFSEH2{i,1},2);
end

bestIndicesSSCFSREH2=bestIndicesSSCFSREH;
for i = 2:6
    bestIndicesSSCFSREH2{i,1}(1,FSRSenDel(1,i))= FSRSenOptimal(1,i); 
    bestIndicesSSCFSREH2{i,1}=sort(bestIndicesSSCFSREH2{i,1},2);
end

optimalSSCFSEH2 =optimalSSCFSEH2';
optimalSSCFSREH2 =optimalSSCFSREH2';

%% 
save OptResSSCExhFS2.mat bestIndicesSSCFSEH2 optimalSSCFSEH2
save OptResSSCExhFSR2.mat bestIndicesSSCFSREH2 optimalSSCFSREH2

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
