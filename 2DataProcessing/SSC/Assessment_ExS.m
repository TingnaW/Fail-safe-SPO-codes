clear;
load OptResSSCExh_Mode1 optimalSSCEH bestIndicesSSCEH
optimalSSC(:,1)= optimalSSCEH;
bestIndicesSSC{1,1}= bestIndicesSSCEH;

load OptResSSCExh_Mode2 optimalSSCEH bestIndicesSSCEH
optimalSSC(:,2)= optimalSSCEH;
bestIndicesSSC{1,2}= bestIndicesSSCEH;

load OptResSSCExh_Mode4 optimalSSCEH bestIndicesSSCEH
optimalSSC(:,3)= optimalSSCEH;
bestIndicesSSC{1,3}= bestIndicesSSCEH;

%% ======= Compare the SSC results corresponding to different modes =======
NumCom=size(optimalSSC,1);
X = 1:1:NumCom;

figure(1)
C = lines(3);FS =15;FS2=15;

p1 = plot(X,optimalSSC(:,1),'-s','color',C(1,:),'LineWidth',1.5);
hold on
p2 = plot(X,optimalSSC(:,2),'-o','color',C(2,:),'LineWidth',1.5);
p3 = plot(X,optimalSSC(:,3),'->','color',C(3,:),'LineWidth',1.5);

xlabel('Number of sensors');
ylabel('SSC')
set(gca,'FontSize',FS,'fontname','times')
ax = gca;
ax.XAxis.FontSize = FS;
ax.YAxis.FontSize = FS;
legend([p1,p2,p3],'Mode 1','Mode 2','Mode 3','FontSize',FS2,'fontname','times','location','northwest')
