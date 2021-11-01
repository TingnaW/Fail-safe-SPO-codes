clear;
FS=15;
NumSen=8; % Set the number of selected sensors

%% ======== Extract the coordinates of sensors ============================
% addpath('add the location of the folder ModeShapeData')  
addpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/ModeShapeData')
load ModeShapeF_T15UD2.mat

% DOF coordinates in space
x= modal_data.Geometry.XYZ(1:36,2);
geom.x(:,1) = modal_data.Geometry.XYZ(1:12,2);
geom.x(:,2) = modal_data.Geometry.XYZ(13:24,2);
geom.x(:,3) = modal_data.Geometry.XYZ(25:36,2);

y =modal_data.Geometry.XYZ(1:36,1);
geom.y(:,1) = modal_data.Geometry.XYZ(1:12,1);
geom.y(:,2) = modal_data.Geometry.XYZ(13:24,1);
geom.y(:,3) = modal_data.Geometry.XYZ(25:36,1);

% rmpath('add the location of the folder ModeShapeData')
rmpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/ModeShapeData')

%% ===== ADPR corresponding to the data set used for the modal identification =======
load Coe4Modes.mat Coe4Modes NF
NFFIM =NF; clear NF;
Coe4ModeM =10^3*Coe4Modes{1,4};

DOFs = size(Coe4ModeM,1);
MSNum = size(NFFIM,2);
ADPR=[];
for i=1:DOFs
    for k=1:MSNum
        ADPR(i,k)= Coe4ModeM(i,k)^2/(NFFIM(1,k)*2*pi);
    end
end
ADPRFIMOri= sum(ADPR,2);
ADPRSSCFIM = normalize(ADPRFIMOri,'range');  

Maxk=[];
[MaxK,~] = maxk(ADPRSSCFIM,NumSen);
ADPRFIMshapes(:,1) = ADPRSSCFIM(1:12,1); 
ADPRFIMshapes(:,2) = ADPRSSCFIM(13:24,1);
ADPRFIMshapes(:,3) = ADPRSSCFIM(25:36,1);

row=[];column=[];
for i =1:NumSen
    [row(i), column(i)]= find(MaxK(i)==ADPRFIMshapes);
end

%
figure(1)
stem3(geom.x,geom.y,ADPRFIMshapes,'bo','MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','b')
axis equal
hold on
for i=1:6
    stem3(geom.x(row(i), column(i)),geom.y(row(i), column(i)),ADPRFIMshapes(row(i), column(i)),'ro','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r')
end

for i=7:8
    stem3(geom.x(row(i), column(i)),geom.y(row(i), column(i)),ADPRFIMshapes(row(i), column(i)),'mo','MarkerSize',5,'MarkerEdgeColor','m','MarkerFaceColor','m')
end

set(gca,'FontSize',FS,'FontName','times')
ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
set(gca,'XTick',[],'YTick',[])
hold off

%% ===== ADPR corresponding to the data set used for the damage detection ===========

load Coe4ModesNorDam.mat X Y NF
NFSSC =NF; clear NF;
MNum = 2; % Select the mode shape, 1:{1,1}; 2:{1,2}; 3{1,3}; 4{1,4};
CoeModeS = 10^3*X{1,MNum}; 

NumObe = size(CoeModeS,1);
DoFs = size(CoeModeS,2);
ADPR=[];
for i=1:DoFs
    for k=1:NumObe
        ADPR(k,i)= CoeModeS(k,i)^2/(NFSSC(k,2)*2*pi); % (NF(k,2); 2nd NF
    end
end
ADPRSSCOri= sum(ADPR,1);
ADPRSSC= normalize(ADPRSSCOri,'range')'; 

Maxk=[];Indall=[];
[MaxK, Indall] = maxk(ADPRSSC,NumSen);
ADPRSSCshapes(:,1) = ADPRSSC(1:12,1); 
ADPRSSCshapes(:,2) = ADPRSSC(13:24,1);
ADPRSSCshapes(:,3) = ADPRSSC(25:36,1);

row=[];column=[];
for i =1:NumSen
    [row(i), column(i)]= find(MaxK(i)==ADPRSSCshapes);
end

figure(2)
stem3(geom.x,geom.y,ADPRSSCshapes,'bo','MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','b')
axis equal
hold on
for i=1:6 
    stem3(geom.x(row(i), column(i)),geom.y(row(i), column(i)),ADPRSSCshapes(row(i), column(i)),'ro','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r')
end

for i=7:8 
    stem3(geom.x(row(i), column(i)),geom.y(row(i), column(i)),ADPRSSCshapes(row(i), column(i)),'mo','MarkerSize',5,'MarkerEdgeColor','m','MarkerFaceColor','m')
end

set(gca,'FontSize',FS,'FontName','times')
ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
set(gca,'XTick',[],'YTick',[])
hold off