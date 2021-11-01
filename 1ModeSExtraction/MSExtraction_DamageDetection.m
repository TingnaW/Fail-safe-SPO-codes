clear;

%% ======== Extract the mode shapes for the next step  ====================
% addpath('add the location of the folder ModeShapeData')  
addpath('/Volumes/GoogleDrive/My Drive/2Data Cluster/Github/Fail-safe/ModeShapeData')

load ModeShapeF_T5UD2.mat
modaldata{1}=modal_data;
load ModeShapeF_T5M1_2.mat
modaldata{2}=modal_data;
load ModeShapeF_T5M2_2.mat
modaldata{3}=modal_data;
load ModeShapeF_T5M3_2.mat
modaldata{4}=modal_data;

load ModeShapeF_T10UD2.mat
modaldata{5}=modal_data;
load ModeShapeF_T10M1_2.mat
modaldata{6}=modal_data;
load ModeShapeF_T10M2_2.mat
modaldata{7}=modal_data;
load ModeShapeF_T10M3_2.mat
modaldata{8}=modal_data;

load ModeShapeF_T15UD2.mat
modaldata{9}=modal_data;
load ModeShapeF_T15M1_2.mat
modaldata{10}=modal_data;
load ModeShapeF_T15M2_2.mat
modaldata{11}=modal_data;
load ModeShapeF_T15M3_2.mat
modaldata{12}=modal_data;

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

for i = 1:MSNum 
    for j = 1:TNum
        X{1,i}(j,:)= Coe4Modes{1,j}(:,i)';
    end
end

y =[0;1;2;3];
y= repmat(y,TNum/4,1); % TNum/4 is equal to the number of the controlled temperature.

Y = dummyvar(categorical(y));
Y = Y(:, 2:end);  

%% ======== Save the required data  =======================================
save('Coe4ModesNorDam.mat', 'X', 'Y','NF')
