clear;
SN = 36;

%% ================ Select 2 sensors ======================================
SelectedSN =2;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndicesFSR=[];optimalSSCIndexFSR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PopulationSize',50,'PlotFcn',@gaplotbestf);
for j = 1:10
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

[x,fval,exitflag] = ga(@SSC_GAFucFSR,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndicesFSR(j,:)=[x];
optimalSSCIndexFSR(j)=abs(fval);    
end

[B,I]= sort(optimalSSCIndexFSR);
bestIndicesSSCFSRGA{SelectedSN,1} = bestIndicesFSR(I(1,end),:)
optimalSSCFSRGA(SelectedSN,1) = B(1,end)

%% ================ Select 3 sensors ======================================
SelectedSN =3;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndicesFSR=[];optimalSSCIndexFSR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PopulationSize',50,'PlotFcn',@gaplotbestf);
for j = 1:10
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

[x,fval,exitflag] = ga(@SSC_GAFucFSR,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndicesFSR(j,:)=[x];
optimalSSCIndexFSR(j)=abs(fval);    
end

[B,I]= sort(optimalSSCIndexFSR);
bestIndicesSSCFSRGA{SelectedSN,1} = bestIndicesFSR(I(1,end),:)
optimalSSCFSRGA(SelectedSN,1) = B(1,end)

%% ================ Select 4 sensors ======================================
SelectedSN =4;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndicesFSR=[];optimalSSCIndexFSR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PopulationSize',50,'PlotFcn',@gaplotbestf);
for j = 1:20
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

[x,fval,exitflag] = ga(@SSC_GAFucFSR,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndicesFSR(j,:)=[x];
optimalSSCIndexFSR(j)=abs(fval);    
end

[B,I]= sort(optimalSSCIndexFSR);
bestIndicesSSCFSRGA{SelectedSN,1} = bestIndicesFSR(I(1,end),:)
optimalSSCFSRGA(SelectedSN,1) = B(1,end)

%% ================ Select 5 sensors ======================================
SelectedSN =5;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndicesFSR=[];optimalSSCIndexFSR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PopulationSize',50,'PlotFcn',@gaplotbestf);
for j = 1:10
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

[x,fval,exitflag] = ga(@SSC_GAFucFSR,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndicesFSR(j,:)=[x];
optimalSSCIndexFSR(j)=abs(fval);    
end

[B,I]= sort(optimalSSCIndexFSR);
bestIndicesSSCFSRGA{SelectedSN,1} = bestIndicesFSR(I(1,end),:)
optimalSSCFSRGA(SelectedSN,1) = B(1,end)

%% ================ Select 6 sensors ======================================
SelectedSN =6;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndicesFSR=[];optimalSSCIndexFSR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PopulationSize',50,'PlotFcn',@gaplotbestf);
for j = 1:20
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

[x,fval,exitflag] = ga(@SSC_GAFucFSR,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndicesFSR(j,:)=[x];
optimalSSCIndexFSR(j)=abs(fval);    
end

[B,I]= sort(optimalSSCIndexFSR);
bestIndicesSSCFSRGA{SelectedSN,1} = bestIndicesFSR(I(1,end),:)
optimalSSCFSRGA(SelectedSN,1) = B(1,end)

%% ================ Select 7 sensors ======================================
SelectedSN =7;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndicesFSR=[];optimalSSCIndexFSR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PopulationSize',50,'PlotFcn',@gaplotbestf);
for j = 1:20
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

[x,fval,exitflag] = ga(@SSC_GAFucFSR,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndicesFSR(j,:)=[x];
optimalSSCIndexFSR(j)=abs(fval);    
end

[B,I]= sort(optimalSSCIndexFSR);
bestIndicesSSCFSRGA{SelectedSN,1} = bestIndicesFSR(I(1,end),:)
optimalSSCFSRGA(SelectedSN,1) = B(1,end)

%% ================ Select 8 sensors ======================================
SelectedSN =8;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndicesFSR=[];optimalSSCIndexFSR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PopulationSize',50,'PlotFcn',@gaplotbestf);
for j = 1:20
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

[x,fval,exitflag] = ga(@SSC_GAFucFSR,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndicesFSR(j,:)=[x];
optimalSSCIndexFSR(j)=abs(fval);    
end

[B,I]= sort(optimalSSCIndexFSR);
bestIndicesSSCFSRGA{SelectedSN,1} = bestIndicesFSR(I(1,end),:)
optimalSSCFSRGA(SelectedSN,1) = B(1,end)

%% ======== Save the results  =============================================
save('OptResSSCGAFSR','optimalSSCFSRGA','bestIndicesSSCFSRGA');

%% ======== Function  =====================================================
function MEd = SSC_GAFucFSR(x)
load Coe4ModesNorDam.mat
CoeModeS = 10^6*X{1,2}; % Select the mode shape, 1:{1,1}; 2:{1,2}; 3{1,3}; 4{1,4};

n=size(x,2);
CoeSelM =[];
    for i=1:n
        CoeSelM = [CoeSelM CoeModeS(:,x(1,i))];  
    end
    
    CoeSelMMid=[];CCC=[];
    for i=1:n
        CoeSelMMid = CoeSelM;
        CoeSelMMid(:,i) = [];  

        [~,~,CCC] = canoncorr(CoeSelMMid,Y);          
        MEdMid(1,i)= sum(CCC.^2);   
    end  
    
[~,Ind] = min(MEdMid);
MEdMid(Ind) = [];
 
MEd= -min(MEdMid); 
end
