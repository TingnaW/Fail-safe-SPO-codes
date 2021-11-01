clear;
load Coe4ModesNorDam.mat X Y NF

MNum = 2; % Select the mode shape, 1:{1,1}; 2:{1,2}; 3{1,3}; 4{1,4};
CoeModeS = 10^6*X{1,MNum}; 

SNTol= 36;
CandiSL = 1:1:SNTol;

%% =========================== Load results ===============================
load OptResSSCADPRExhFS.mat bestIndicesSSCADPRFSEH optimalSSCADPRFSEH
load OptResSSCADPRExhFSR.mat bestIndicesSSCADPRFSREH optimalSSCADPRFSREH

%% *************************** Two Sensor locations ***********************
NumSen=2;

% GAR 
C =[];C = bestIndicesSSCADPRFSEH{NumSen,1};
[optimalSSCADPRFSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCADPRMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y,NF);

% GAFSR 
C =[];C = bestIndicesSSCADPRFSREH{NumSen,1};
[optimalSSCADPRFSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCADPRMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y,NF);

%% *************************** Three Sensor locations *********************
NumSen=3;

% GAR 
C =[];C = bestIndicesSSCADPRFSEH{NumSen,1};
[optimalSSCADPRFSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCADPRMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y,NF);

% GAFSR 
C =[];C = bestIndicesSSCADPRFSREH{NumSen,1};
[optimalSSCADPRFSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCADPRMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y,NF);

%% *************************** Four Sensor locations **********************
NumSen=4;

% GAR 
C =[];C = bestIndicesSSCADPRFSEH{NumSen,1};
[optimalSSCADPRFSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCADPRMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y,NF);

% GAFSR 
C =[];C = bestIndicesSSCADPRFSREH{NumSen,1};
[optimalSSCADPRFSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCADPRMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y,NF);

%% *************************** Five Sensor locations **********************
NumSen=5;

% GAR 
C =[];C = bestIndicesSSCADPRFSEH{NumSen,1};
[optimalSSCADPRFSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCADPRMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y,NF);

% GAFSR 
C =[];C = bestIndicesSSCADPRFSREH{NumSen,1};
[optimalSSCADPRFSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCADPRMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y,NF);

%% *************************** Six Sensor locations ***********************
NumSen=6;

% GAR 
C =[];C = bestIndicesSSCADPRFSEH{NumSen,1};
[optimalSSCADPRFSEH2(1,NumSen), FSSenDel(1,NumSen),FSSenOptimal(1,NumSen),FSSenAllPossi{1,NumSen}]= SSCADPRMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y,NF);

% GAFSR 
C =[];C = bestIndicesSSCADPRFSREH{NumSen,1};
[optimalSSCADPRFSREH2(1,NumSen), FSRSenDel(1,NumSen),FSRSenOptimal(1,NumSen),FSRSenAllPossi{1,NumSen}]= SSCADPRMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y,NF);

%% ======== Save the results  =============================================
bestIndicesSSCADPRFSEH2=bestIndicesSSCADPRFSEH;
for i = 2:6
    bestIndicesSSCADPRFSEH2{i,1}(1,FSSenDel(1,i))= FSSenOptimal(1,i); 
    bestIndicesSSCADPRFSEH2{i,1}=sort(bestIndicesSSCADPRFSEH2{i,1},2);
end

bestIndicesSSCADPRFSREH2=bestIndicesSSCADPRFSREH;
for i = 2:6
    bestIndicesSSCADPRFSREH2{i,1}(1,FSRSenDel(1,i))= FSRSenOptimal(1,i); 
    bestIndicesSSCADPRFSREH2{i,1}=sort(bestIndicesSSCADPRFSREH2{i,1},2);
end

optimalSSCADPRFSEH2 =optimalSSCADPRFSEH2';
optimalSSCADPRFSREH2 =optimalSSCADPRFSREH2';

%%
save OptResSSCADPRExhFS2.mat bestIndicesSSCADPRFSEH2 optimalSSCADPRFSEH2
save OptResSSCADPRExhFSR2.mat bestIndicesSSCADPRFSREH2 optimalSSCADPRFSREH2

%% ======== Function 1  ===================================================
function [SSCADPRMaxFS,I,OptimalS,CandiSL0]= SSCADPRMax_FucFS(NumSen,C,CoeModeS,CandiSL,Y,NF)
    NumObe = size(CoeModeS,1);
    DoFs = size(CoeModeS,2);
    for i=1:DoFs
        for k=1:NumObe
            ADPR(k,i)= CoeModeS(k,i)^2/(NF(k,2)*2*pi); % (NF(k,2); 2nd NF
        end
    end
    ADPRDOF= sum(ADPR,1);
    ADPRDOF = normalize(ADPRDOF,'range'); 

    CoeSelM =[]; ADPRM=[];
    for i=1:NumSen
        CoeSelM = [CoeSelM CoeModeS(:,C(1,i))];  
        ADPRM = [ADPRM ADPRDOF(1,C(1,i))];
    end

    for i=1:NumSen
        CoeSelMMid = CoeSelM;
        CoeSelMMid(:,i) = [];
        ADPRMMid = ADPRM;
        ADPRMMid(:,i) = []; 
        
        [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
        MEdMid(1,i)= sum(CCC.^2)*sum(ADPRMMid);   
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
        CoeSelM =[]; ADPRM=[];
        for j=1:NumSen
            CoeSelM = [CoeSelM CoeModeS(:,C(1,j))];  
            ADPRM = [ADPRM ADPRDOF(1,C(1,j))];
        end
        
        for j=1:NumSen
            CoeSelMMid = CoeSelM;
            CoeSelMMid(:,j) = [];
            ADPRMMid = ADPRM;
            ADPRMMid(:,j) = []; 
            
            [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
            MEdMid(1,j)= sum(CCC.^2)*sum(ADPRMMid); 
        end  
        
        if min(MEdMid)<SSCFS
           CandiSL0(:,CandiSL0==CandiSL(1,i))=[];
        else
           [~,~,CCC] = canoncorr(CoeSelM,Y);          
           SSC(1,i)=sum(CCC.^2)*sum(ADPRM);  
        end
    end 
    [SSCADPRMaxFS,Ind] = max(SSC);
    OptimalS = CandiSL(1,Ind);
end

%% ======== Function 2  ===================================================
function [SSCADPRMaxFSR,Idel,OptimalS,CandiSL0]= SSCADPRMax_FucFSR(NumSen,C,CoeModeS,CandiSL,Y,NF)
    NumObe = size(CoeModeS,1);
    DoFs = size(CoeModeS,2);
    for i=1:DoFs
        for k=1:NumObe
            ADPR(k,i)= CoeModeS(k,i)^2/(NF(k,2)*2*pi); % (NF(k,2); 2nd NF
        end
    end
    ADPRDOF= sum(ADPR,1);
    ADPRDOF = normalize(ADPRDOF,'range'); 

    CoeSelM =[]; ADPRM=[];
    for i=1:NumSen
        CoeSelM = [CoeSelM CoeModeS(:,C(1,i))];  
        ADPRM = [ADPRM ADPRDOF(1,C(1,i))];
    end

    for i=1:NumSen
        CoeSelMMid = CoeSelM;
        CoeSelMMid(:,i) = [];
        ADPRMMid = ADPRM;
        ADPRMMid(:,i) = [];         
        
        [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
        MEdMid(1,i)= sum(CCC.^2)*sum(ADPRMMid);  
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
        CoeSelM =[]; ADPRM=[];
        for j=1:NumSen
            CoeSelM = [CoeSelM CoeModeS(:,C(1,j))];  
            ADPRM = [ADPRM ADPRDOF(1,C(1,j))];
        end
        
        for j=1:NumSen
            CoeSelMMid = CoeSelM;
            CoeSelMMid(:,j) = [];
            ADPRMMid = ADPRM;
            ADPRMMid(:,j) = [];       
            
            [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
            MEdMid(1,j)= sum(CCC.^2)*sum(ADPRMMid); 
        end  
        
        [~,IRep] = min(MEdMid);
        MEdMid(:,IRep) = [];   
        SSCFSR2 = min(MEdMid);
        if SSCFSR2<SSCFSR
           CandiSL0(:,CandiSL0==CandiSL(1,i))=[];
        else
           [~,~,CCC] = canoncorr(CoeSelM,Y);          
           SSC(1,i)=sum(CCC.^2)*sum(ADPRM);    
        end
    end 
    [SSCADPRMaxFSR,Ind] = max(SSC);
    OptimalS = CandiSL(1,Ind);
end
