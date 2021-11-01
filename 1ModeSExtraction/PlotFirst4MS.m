clear;
FS= 13; % fontsize

%% ======== Extract the mode shapes for the health state ==================
% addpath('add the location of the folder ModeShapeData')  
addpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/ModeShapeData')

load ModeShapeF_T0UD2.mat
modaldata{1}=modal_data;
load ModeShapeF_T5UD2.mat
modaldata{2}=modal_data;
load ModeShapeF_T10UD2.mat
modaldata{3}=modal_data;
load ModeShapeF_T15UD2.mat
modaldata{4}=modal_data;
load ModeShapeF_T20UD2.mat
modaldata{5}=modal_data;
load ModeShapeF_T25UD2.mat
modaldata{6}=modal_data;

% rmpath('add the location of the folder ModeShapeData')
rmpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/ModeShapeData')

NFCol=[1 2 3 4]; % Select mode shapes
MSNum=size(NFCol,2); % The number of mode shapes
TNum = size(modaldata,2);

NF =zeros(TNum,MSNum); Coe4Modes =[];
for j = 1:TNum 
    for i = 1:MSNum
        NF(j,:)=modaldata{j}.modeFreq(1,NFCol);
        Coe4Modes{1,j}(:,i)=modaldata{j}.modeShape{1,NFCol(1,i)}(:,3); 
    end
end

%% ======== Plot the mode shapes for the health state =====================
DOFs=size(Coe4Modes{1,1},1);  % The degree of freedoms

% DOF coordinates in space
geom.x{1}(:,1) = modaldata{1}.Geometry.XYZ(1:12,2);
geom.x{1}(:,2) = modaldata{1}.Geometry.XYZ(13:24,2);
geom.x{1}(:,3) = modaldata{1}.Geometry.XYZ(25:36,2);

geom.y{1}(:,1) = modaldata{1}.Geometry.XYZ(1:12,1);
geom.y{1}(:,2) = modaldata{1}.Geometry.XYZ(13:24,1);
geom.y{1}(:,3) = modaldata{1}.Geometry.XYZ(25:36,1);

% reorder mode shapes according to DOFs
for j =1:TNum
    for i =1:MSNum
        modeshapes{1,j}.sub{i,1}(:,1) = Coe4Modes{1,j}(1:12,i); % 0.5*i*10^4
        modeshapes{1,j}.sub{i,1}(:,2) = Coe4Modes{1,j}(13:24,i);
        modeshapes{1,j}.sub{i,1}(:,3) = Coe4Modes{1,j}(25:36,i);
    end
end

for j =1:TNum
    figure(j)
    for i = 1:MSNum
        titles = ['Mode ',num2str(i), ', freq = ',num2str(NF(j,i)), ' Hz'];
        subplot(2,2,i)
        surf(geom.x{1,1}, geom.y{1,1}, modeshapes{1,j}.sub{i,1})
        colormap autumn
        title(titles,'fontname','times','fontsize',FS)
        zticks(0)
        set(gca,'fontSize',FS,'fontname','times') 
    end
end