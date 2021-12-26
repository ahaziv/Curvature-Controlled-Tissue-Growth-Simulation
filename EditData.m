% this code comes to edit .mat results in case of some minor alterations
% are required.
clear all; close all; clc;
addpath('D:\Documents\MATLAB\Thesis\3D CCTG\Output data')

ShapeType = 'Scaffold_27_Units';        % NormRatio = 3372.53121666825
% ShapeType = 'Scaffold_Unit_Plus_Plus';        % NormRatio = 8090.54723387182

%% loading
Loadfilename = strcat(ShapeType,'.mat');
ShapeData = importdata(Loadfilename);

%% editing
NormRatio = 3372.53121666825;
Lambda = ShapeData{4};
Lambda = Lambda/(NormRatio^2);
ShapeData{4} = Lambda;

%% saving
Savefilename = strcat('D:\Documents\MATLAB\Thesis\3D CCTG\Output data\',ShapeType);
file = ShapeData;
save(Savefilename, 'file')