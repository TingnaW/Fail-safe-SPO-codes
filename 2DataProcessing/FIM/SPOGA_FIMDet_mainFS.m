clear;
SN = 36;

%% ============== Select 4 optimal sensors ================================
NumSen = 4;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndicesR=[];optimalFIMIndexR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PlotFcn',@gaplotbestf);
for j = 1:10
    nVar = NumSen; % Here, you can set the number of selected sensors  
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

[x,fval,exitflag] = ga(@FIM_GAFucFS,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndicesR(j,:)=[x];
optimalFIMIndexR(j)=abs(fval);       
end

[B,I]= sort(optimalFIMIndexR);
bestIndicesFIMGAFS{NumSen,1} = bestIndicesR(I(1,end),:)
optimalFIMGAFS(NumSen,1) = B(1,end)

%% ============== Select 5 optimal sensors ================================
NumSen = 5;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndicesR=[];optimalFIMIndexR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PlotFcn',@gaplotbestf);
for j = 1:10
    nVar = NumSen; % Here, you can set the number of selected sensors  
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

[x,fval,exitflag] = ga(@FIM_GAFucFS,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndicesR(j,:)=[x];
optimalFIMIndexR(j)=abs(fval);       
end

[B,I]= sort(optimalFIMIndexR);
bestIndicesFIMGAFS{NumSen,1} = bestIndicesR(I(1,end),:)
optimalFIMGAFS(NumSen,1) = B(1,end)

%% ============== Select 6 optimal sensors ================================
NumSen = 6;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndicesR=[];optimalFIMIndexR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PlotFcn',@gaplotbestf);
for j = 1:10
    nVar = NumSen; % Here, you can set the number of selected sensors  
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

[x,fval,exitflag] = ga(@FIM_GAFucFS,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndicesR(j,:)=[x];
optimalFIMIndexR(j)=abs(fval);       
end

[B,I]= sort(optimalFIMIndexR);
bestIndicesFIMGAFS{NumSen,1} = bestIndicesR(I(1,end),:)
optimalFIMGAFS(NumSen,1) = B(1,end)

%% ============== Select 7 optimal sensors ================================
NumSen = 7;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndicesR=[];optimalFIMIndexR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PlotFcn',@gaplotbestf);
for j = 1:10
    nVar = NumSen; % Here, you can set the number of selected sensors  
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

[x,fval,exitflag] = ga(@FIM_GAFucFS,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndicesR(j,:)=[x];
optimalFIMIndexR(j)=abs(fval);       
end

[B,I]= sort(optimalFIMIndexR);
bestIndicesFIMGAFS{NumSen,1} = bestIndicesR(I(1,end),:)
optimalFIMGAFS(NumSen,1) = B(1,end)

%% ============== Select 8 optimal sensors ================================
NumSen = 8;

lb =[]; ub=[]; IntCon =[]; C =[]; D=[];x =[]; fval =[];
bestIndicesR=[];optimalFIMIndexR=[]; B=[]; I=[];
opts = optimoptions('ga','ConstraintTolerance',1e-30,'FunctionTolerance',1e-30,'PlotFcn',@gaplotbestf);
for j = 1:10
    nVar = NumSen; % Here, you can set the number of selected sensors  
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

[x,fval,exitflag] = ga(@FIM_GAFucFS,nVar,C,D,[],[],...
                       lb,ub,[],IntCon,opts);

bestIndicesR(j,:)=[x];
optimalFIMIndexR(j)=abs(fval);       
end

[B,I]= sort(optimalFIMIndexR);
bestIndicesFIMGAFS{NumSen,1} = bestIndicesR(I(1,end),:)
optimalFIMGAFS(NumSen,1) = B(1,end)

%% ======== Save the results  =============================================
save('OptResFIMGAFS','optimalFIMFSGA','bestIndicesFIMFSGA');
    
%% ======== Function  =====================================================
function MEd = FIM_GAFucFS(x)
    load Coe4Modes.mat
    Coe4ModeM = 10^3*Coe4Modes{1,4}; 
    n=size(x,2);
    CoeSecM =[];
    for i=1:n
        CoeSecM = [CoeSecM;Coe4ModeM(x(1,i),:)];  
    end
    
    CoeSecMMid =[];
    for i=1:n
        CoeSecMMid = CoeSecM;
        CoeSecMMid(i,:) = [];
        MEdMid(1,i)= det(CoeSecMMid.'*CoeSecMMid);
    end    
MEd = -min(MEdMid);
end

