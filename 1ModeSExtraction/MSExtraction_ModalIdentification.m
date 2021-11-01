clear;

%% ======== Extract the mode shapes for the next step  ====================
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

NFCol=[1 2 4]; % Select the mode shapes going to be used.
MSNum=size(NFCol,2); % The number of mode shapes.
TNum = size(modaldata,2); % The number of the controlled temperature.

% Reorganise the data to put them into a matrix
NF =zeros(TNum,MSNum); Coe4Modes =[];
for j = 1:TNum 
    for i = 1:MSNum
        NF(j,:)=modaldata{j}.modeFreq(1,NFCol);
        Coe4Modes{1,j}(:,i)=modaldata{j}.modeShape{1,NFCol(1,i)}(:,3); 
    end
end

%% ======== Save the required data  =======================================
save Coe4Modes.mat Coe4Modes NF
