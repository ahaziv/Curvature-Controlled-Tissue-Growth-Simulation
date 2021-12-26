% the fast plotter 
clear all; close all; clc;
addpath('D:\Documents\MATLAB\Thesis\3D CCTG\Output data')

ShapeType = 'Banana';
% ShapeType = 'Scaffold_27_Units';
PlotRate = 0.25;
FlagPartView = 0;
FlagDispInitGeom = 0;
FlagBoundary = 0;

filename = strcat(ShapeType,'.mat');
ShapeData = importdata(filename);
Faces = ShapeData{1}; Vertices = ShapeData{2}; ColorVec = ShapeData{3}; Lambda = ShapeData{4}; VGrowth = ShapeData{5}; Vtot = ShapeData{6};

if FlagDispInitGeom
    InitGeom = {Faces{1} Vertices{1}};
else
    InitGeom = [];
end

if FlagPartView%length(Faces{1}(1,:))>3
%     FlagPartView = 1;
    inbVertices = unique(reshape(Faces{1}(Faces{1}(:,4)==1,1:3),[numel(Faces{1}(Faces{1}(:,4)==1,1:3)),1]));
%     UnitVol = (max(Vertices{1}(inbVertices,1))-min(Vertices{1}(inbVertices,1)))*(max(Vertices{1}(inbVertices,2))-min(Vertices{1}(inbVertices,2)))*(max(Vertices{1}(inbVertices,3))-min(Vertices{1}(inbVertices,3))); 
    maxis = 1.3*max(max(abs(Vertices{1}(inbVertices,1:3))));
else
%     UnitVol = (max(Vertices{1}(:,1))-min(Vertices{1}(:,1)))*(max(Vertices{1}(:,2))-min(Vertices{1}(:,2)))*(max(Vertices{1}(:,3))-min(Vertices{1}(:,3)));
    maxis = 1.3*max(max(abs(Vertices{1}(:,1:3))));
end
minmaxColor = [min(ColorVec{1}) max(abs(ColorVec{1}))];
TimeVec = zeros(length(Lambda),1);
TimeVec(1) = Lambda(1);
for jj = 2:length(Lambda)
    TimeVec(jj) = TimeVec(jj-1)+Lambda(jj);
end
%% open a figure
fig_h = figure('name',ShapeType,'numbertitle','off','color',[0.75 0.75 0.75]);
ax = axes('DataAspectRatio', [1,1,1]);
set(gcf, 'OuterPosition', get(0, 'Screensize'));
axis tight
%% plotting
critTimePoints = [1 25 45 45];

for ii = 1:length(Faces)
    ColorVec{ii}(:) = -2;
    PlotMesh(Faces{ii}, Vertices{ii}, ColorVec{ii}, maxis, minmaxColor, FlagBoundary, FlagPartView, InitGeom) 
    pause(PlotRate)
    if ismember(ii, critTimePoints)
        abc=1;
    end
end

%% plotting the volume fill
figure
NormVGrowthRatio = VGrowth./(Lambda*Vtot);
plot(TimeVec,NormVGrowthRatio)
title('volume ratio filling rate')
xlabel('Time')
ylabel('V_f_i_l_l/(V_t_o_t*\Deltat)')

figure
NormVGrowthRatio = VGrowth./(Lambda*Vtot);
plot(TimeVec,NormVGrowthRatio)
title('volume ratio filling rate')
xlabel('Time')
ylabel('V_f_i_l_l/(V_t_o_t*\Deltat)')

text = strcat('Total Volume growth Ratio: ',num2str(sum(VGrowth)/Vtot));
disp(text)
