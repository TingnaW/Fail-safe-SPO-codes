clear;
SN = 36;

%% ================ Select 2 sensors ======================================
SelectedSN =2;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndices=[];optimalSSCADPRIndexFSR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PopulationSize',50,'PlotFcn',@gaplotbestf);
for j = 1:5
    nVar = SelectedSN; % Here, you can set the number of selected sensors  
    fprintf(['Select ',num2str(nVar),' optimal sensors...\n']);
    lb = ones(1,nVar);
    ub = SN*ones(1,nVar);
    IntCon = [1:nVar];

    nVar2 = nVar-1;
    C = zeros(nVar2,nVar);
    D = -ones(nVar2,1);

    for k =1 : nVar2
        C(k,k)=1;
        C(k,k+1)=-1;
    end

[x,fval,exitflag] = ga(@SSCADPR_GAFucFSR,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndices(j,:)=[x];
optimalSSCADPRIndexFSR(j)=abs(fval);    
end

[B,I]= sort(optimalSSCADPRIndexFSR);
bestIndicesSSCADPRFSRGA{SelectedSN,1} = bestIndices(I(1,end),:)
optimalSSCADPRFSRGA(SelectedSN,1) = B(1,end)

%% ================ Select 3 sensors ======================================
SelectedSN =3;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndices=[];optimalSSCADPRIndexFSR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PopulationSize',50,'PlotFcn',@gaplotbestf);
for j = 1:10
    nVar = SelectedSN;  
    fprintf(['Select ',num2str(nVar),' optimal sensors...\n']);
    lb = ones(1,nVar);
    ub = SN*ones(1,nVar);
    IntCon = [1:nVar];

    nVar2 = nVar-1;
    C = zeros(nVar2,nVar);
    D = -ones(nVar2,1);

    for k =1 : nVar2
        C(k,k)=1;
        C(k,k+1)=-1;
    end

[x,fval,exitflag] = ga(@SSCADPR_GAFucFSR,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndices(j,:)=[x];
optimalSSCADPRIndexFSR(j)=abs(fval);    
end

[B,I]= sort(optimalSSCADPRIndexFSR);
bestIndicesSSCADPRFSRGA{SelectedSN,1} = bestIndices(I(1,end),:)
optimalSSCADPRFSRGA(SelectedSN,1) = B(1,end)

%% ================ Select 4 sensors ======================================
SelectedSN =4;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndices=[];optimalSSCADPRIndexFSR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PopulationSize',50,'PlotFcn',@gaplotbestf);
for j = 1:10
    nVar = SelectedSN;  
    fprintf(['Select ',num2str(nVar),' optimal sensors...\n']);
    lb = ones(1,nVar);
    ub = SN*ones(1,nVar);
    IntCon = [1:nVar];

    nVar2 = nVar-1;
    C = zeros(nVar2,nVar);
    D = -ones(nVar2,1);

    for k =1 : nVar2
        C(k,k)=1;
        C(k,k+1)=-1;
    end

[x,fval,exitflag] = ga(@SSCADPR_GAFucFSR,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndices(j,:)=[x];
optimalSSCADPRIndexFSR(j)=abs(fval);    
end

[B,I]= sort(optimalSSCADPRIndexFSR);
bestIndicesSSCADPRFSRGA{SelectedSN,1} = bestIndices(I(1,end),:)
optimalSSCADPRFSRGA(SelectedSN,1) = B(1,end)

%% ================ Select 5 sensors ======================================
SelectedSN =5;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndices=[];optimalSSCADPRIndexFSR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PopulationSize',50,'PlotFcn',@gaplotbestf);
for j = 1:10
    nVar = SelectedSN;   
    fprintf(['Select ',num2str(nVar),' optimal sensors...\n']);
    lb = ones(1,nVar);
    ub = SN*ones(1,nVar);
    IntCon = [1:nVar];

    nVar2 = nVar-1;
    C = zeros(nVar2,nVar);
    D = -ones(nVar2,1);

    for k =1 : nVar2
        C(k,k)=1;
        C(k,k+1)=-1;
    end

[x,fval,exitflag] = ga(@SSCADPR_GAFucFSR,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndices(j,:)=[x];
optimalSSCADPRIndexFSR(j)=abs(fval);    
end

[B,I]= sort(optimalSSCADPRIndexFSR);
bestIndicesSSCADPRFSRGA{SelectedSN,1} = bestIndices(I(1,end),:)
optimalSSCADPRFSRGA(SelectedSN,1) = B(1,end)

%% ================ Select 6 sensors ======================================
SelectedSN =6;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndices=[];optimalSSCADPRIndexFSR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PopulationSize',50,'PlotFcn',@gaplotbestf);
for j = 1:10
    nVar = SelectedSN;  
    fprintf(['Select ',num2str(nVar),' optimal sensors...\n']);
    lb = ones(1,nVar);
    ub = SN*ones(1,nVar);
    IntCon = [1:nVar];

    nVar2 = nVar-1;
    C = zeros(nVar2,nVar);
    D = -ones(nVar2,1);

    for k =1 : nVar2
        C(k,k)=1;
        C(k,k+1)=-1;
    end

[x,fval,exitflag] = ga(@SSCADPR_GAFucFSR,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndices(j,:)=[x];
optimalSSCADPRIndexFSR(j)=abs(fval);    
end

[B,I]= sort(optimalSSCADPRIndexFSR);
bestIndicesSSCADPRFSRGA{SelectedSN,1} = bestIndices(I(1,end),:)
optimalSSCADPRFSRGA(SelectedSN,1) = B(1,end)

%% ================ Select 7 sensors ======================================
SelectedSN =7;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndices=[];optimalSSCADPRIndexFSR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PopulationSize',50,'PlotFcn',@gaplotbestf);
for j = 1:10
    nVar = SelectedSN;   
    fprintf(['Select ',num2str(nVar),' optimal sensors...\n']);
    lb = ones(1,nVar);
    ub = SN*ones(1,nVar);
    IntCon = [1:nVar];

    nVar2 = nVar-1;
    C = zeros(nVar2,nVar);
    D = -ones(nVar2,1);

    for k =1 : nVar2
        C(k,k)=1;
        C(k,k+1)=-1;
    end

[x,fval,exitflag] = ga(@SSCADPR_GAFucFSR,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndices(j,:)=[x];
optimalSSCADPRIndexFSR(j)=abs(fval);    
end

[B,I]= sort(optimalSSCADPRIndexFSR);
bestIndicesSSCADPRFSRGA{SelectedSN,1} = bestIndices(I(1,end),:)
optimalSSCADPRFSRGA(SelectedSN,1) = B(1,end)

%% ================ Select 8 sensors ======================================
SelectedSN =8;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndices=[];optimalSSCADPRIndexFSR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PopulationSize',50,'PlotFcn',@gaplotbestf);
for j = 1:10
    nVar = SelectedSN;  
    fprintf(['Select ',num2str(nVar),' optimal sensors...\n']);
    lb = ones(1,nVar);
    ub = SN*ones(1,nVar);
    IntCon = [1:nVar];

    nVar2 = nVar-1;
    C = zeros(nVar2,nVar);
    D = -ones(nVar2,1);

    for k =1 : nVar2
        C(k,k)=1;
        C(k,k+1)=-1;
    end

[x,fval,exitflag] = ga(@SSCADPR_GAFucFSR,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndices(j,:)=[x];
optimalSSCADPRIndexFSR(j)=abs(fval);    
end

[B,I]= sort(optimalSSCADPRIndexFSR);
bestIndicesSSCADPRFSRGA{SelectedSN,1} = bestIndices(I(1,end),:)
optimalSSCADPRFSRGA(SelectedSN,1) = B(1,end)

%% ======== Save the results  =============================================
save('OptResSSCADPRGAFSR','optimalSSCADPRFSRGA','bestIndicesSSCADPRFSRGA');

%% ======== Function  ===================================================== 
function MEd = SSCADPR_GAFucFSR(x)
load Coe4ModesNorDam.mat X Y NF
CoeModeS = 10^6*X{1,2}; % Select the mode shape, 1:{1,1}; 2:{1,2}; 3{1,3}; 4{1,4};

NumObe = size(CoeModeS,1);
DoFs = size(CoeModeS,2);
    for i=1:DoFs
        for k=1:NumObe
            ADPR(k,i)= CoeModeS(k,i)^2/(NF(k,2)*2*pi); % (NF(k,2); 2nd NF
        end
    end
    ADPRDOF= sum(ADPR,1);
    ADPRDOF = normalize(ADPRDOF,'range');        

    n=size(x,2);
    CoeSelM =[];ADPRM=[];                        
    for i=1:n
        CoeSelM = [CoeSelM CoeModeS(:,x(1,i))]; 
        ADPRM = [ADPRM ADPRDOF(1,x(1,i))];
    end 
    
    for i=1:n
        CoeSelMMid = CoeSelM;
        CoeSelMMid(:,i) = [];  
        ADPRMMid = ADPRM;
        ADPRMMid(:,i) = []; 

        [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
        MEdMid(1,i)= sum(CCC.^2)*sum(ADPRMMid);     
    end  
[~,Ind] = min(MEdMid);
MEdMid(Ind) = [];
MEd= -min(MEdMid); 
end

