%% Bar chart
% T = [18; 183; 53; 486; 17; 6; 41/3; 70/3; 253/3];
T = [1.77; 0.7292; 0.7274; 0.6238;.1239];

T = T';
% T = sort(T);
figure,
%subplot(1,2,1), 
b = bar(1:size(T),T,'FaceColor','flat');
set(gca,'xtick',[])
set(gca,'xticklabel',[])

set(b(1), 'FaceColor','[0.7 0.8 0.5]');
set(b(2), 'FaceColor','[0.2 0.6 0.44]');
set(b(3), 'FaceColor','[0.75 0.75 0.75]');
set(b(4), 'FaceColor','[.2 .2 .1');
set(b(5), 'FaceColor','[0.7 0.8 0.9]');

lgnd = legend ('tSVD (MATLAB)','WNNM (MATLAB)','Proposed (MATLAB)','SAIST (MATLAB)','BM3D (C++)',...
    'orientation','horizontal','location','northeast','NumColumns',1);
 set(lgnd,'color','none');
 set(lgnd,'edgecolor','none');
set(gca,'FontSize',14);
ylabel('Average CPU time [s]')
grid on;

%% Pie chart

X = [0.9218 5.1438 0.1276 ];
figure,
%subplot(1,2,2),
pie(X)
labels = {'Linesearch','SVDs','Image Update/Other',};
legend(labels)
set(gca,'FontSize',12);
