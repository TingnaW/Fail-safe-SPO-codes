clear;
SN = 36;

%% ================ Select 3 sensors ======================================
SelectedSN = 3;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndices=[];optimalFIMADRPIndex=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PlotFcn',@gaplotbestf);
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

[x,fval,exitflag] = ga(@FIMADPR_GAFuc,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndices(j,:)=[x];
optimalFIMADRPIndex(j)=abs(fval);       
end

[B,I]= sort(optimalFIMADRPIndex);
bestIndicesFIMADPRGA{SelectedSN,1} = bestIndices(I(1,end),:)
optimalFIMADPRGA(SelectedSN,1) = B(1,end)

%% ================ Select 4 sensors ======================================
SelectedSN = 4;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndices=[];optimalFIMADRPIndex=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PlotFcn',@gaplotbestf);
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

[x,fval,exitflag] = ga(@FIMADPR_GAFuc,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndices(j,:)=[x];
optimalFIMADRPIndex(j)=abs(fval);       
end

[B,I]= sort(optimalFIMADRPIndex);
bestIndicesFIMADPRGA{SelectedSN,1} = bestIndices(I(1,end),:)
optimalFIMADPRGA(SelectedSN,1) = B(1,end)

%% ================ Select 5 sensors ======================================
SelectedSN = 5;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndices=[];optimalFIMADRPIndex=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PlotFcn',@gaplotbestf);
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

[x,fval,exitflag] = ga(@FIMADPR_GAFuc,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndices(j,:)=[x];
optimalFIMADRPIndex(j)=abs(fval);       
end

[B,I]= sort(optimalFIMADRPIndex);
bestIndicesFIMADPRGA{SelectedSN,1} = bestIndices(I(1,end),:)
optimalFIMADPRGA(SelectedSN,1) = B(1,end)

%% ================ Select 6 sensors ======================================
SelectedSN = 6;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndices=[];optimalFIMADRPIndex=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PlotFcn',@gaplotbestf);
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

[x,fval,exitflag] = ga(@FIMADPR_GAFuc,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndices(j,:)=[x];
optimalFIMADRPIndex(j)=abs(fval);       
end

[B,I]= sort(optimalFIMADRPIndex);
bestIndicesFIMADPRGA{SelectedSN,1} = bestIndices(I(1,end),:)
optimalFIMADPRGA(SelectedSN,1) = B(1,end)

%% ================ Select 7 sensors ======================================
SelectedSN = 7;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndices=[];optimalFIMADRPIndex=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PlotFcn',@gaplotbestf);
for j = 1:10
    nVar = SelectedSN; % Here, you can set the number of selected sensors  
    fprintf(['Select ',num2str(nVar),' two optimal sensors...\n']);
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

[x,fval,exitflag] = ga(@FIMADPR_GAFuc,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndices(j,:)=[x];
optimalFIMADRPIndex(j)=abs(fval);       
end

[B,I]= sort(optimalFIMADRPIndex);
bestIndicesFIMADPRGA{SelectedSN,1} = bestIndices(I(1,end),:)
optimalFIMADPRGA(SelectedSN,1) = B(1,end)

%% ================ Select 8 sensors ======================================
SelectedSN = 8;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndices=[];optimalFIMADRPIndex=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PlotFcn',@gaplotbestf);
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

[x,fval,exitflag] = ga(@FIMADPR_GAFuc,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndices(j,:)=[x];
optimalFIMADRPIndex(j)=abs(fval);       
end

[B,I]= sort(optimalFIMADRPIndex);
bestIndicesFIMADPRGA{SelectedSN,1} = bestIndices(I(1,end),:)
optimalFIMADPRGA(SelectedSN,1) = B(1,end)

%% ======== Save the results  =============================================
save('OptResFIMADPRGA','optimalFIMADPRGA','bestIndicesFIMADPRGA');
    
%% ======== Function  =====================================================
function MEd = FIMADPR_GAFuc(x)
    load Coe4Modes.mat
    Coe4ModeM = 10^3*Coe4Modes{1,4}; NFS = NF(4,:);
    
    DOFs = size(Coe4ModeM,1);
    MSNum = size(NFS,2);
    for i=1:DOFs
        for k=1:MSNum
            ADPR(i,k)= Coe4ModeM(i,k)^2/(NF(1,k)*2*pi);
        end
    end
    ADPRDOF= sum(ADPR,2);
    ADPRDOF = normalize(ADPRDOF,'range');        
    
    n=size(x,2);
    CoeSecM =[]; ADPRM = [];
    for i=1:n
        CoeSecM = [CoeSecM;Coe4ModeM(x(1,i),:)]; 
        ADPRM = [ADPRM;ADPRDOF(x(1,i),1)];
    end   
MEd = -det(CoeSecM.'*CoeSecM)*sum(ADPRM);
end

