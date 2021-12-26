%% This is the main code for simulating CCTG 
% This code runs the CCTG simulation and will plot each iteration geometry while doing so. In order to work file path must be re-written in the first two lines, so the program may find all relevant functions and data.
clear all; close all; clc;
addpath('C:\Documents\MATLAB\Thesis\3D CCTG\Shapes data')       % file paths
addpath('C:\Documents\MATLAB\Thesis\3D CCTG\pics')              % file paths

%% general parameters
NumIter = 5000;                                                 % number of iteration steps
% TimeLimit = 4.07e-09;
TimeLimit = 1e6;                    % simulation time limit, note that it can greatly vary (orders of magnitude) between geometries

VGrowthRateLim = 0.5*10^-2;         % a growth limit in to detect shape convergence, can be set to zero but then may require a manual termination
LFactor = 0.3;                      % Lambda multiplication factor
DeltaTheta = 0.13;                  % the square of the minimal angle for remeshing 
FlagBoundary = 1;                   % this flag tells the program wether to impose boundary conditions
TimeFactor = 1;%14*10^-14;          % this number represents the biological growth factor [m^2/sec]

FlagAnaSol = 0;                     % can be set to 1 in cases an analytical solution detection is desired
AnaSolution = 0;                    % initial value for variable, should always be set as zero
if FlagAnaSol
    AnalyticParam = struct('GenEps', 0.5, 'SphereEps', 0.2, 'CylinderEps', 0.5);
end

FlagDispInitGeom = 1;               % displays a half transparent plot of the initial geometry througout the run 
FlagVideo = 0;                      % this flag tells the program wether to record the run adn produce a video
if FlagVideo
    TotVideoLength = 20;
    images = cell(NumIter,1);
end
FlagPartView = 0;
FlagSave = 1;
FlagDegradation = 0;

Faces = cell(NumIter,1);
Vertices = cell(NumIter,1);
ShapeType = 'Big_square_scaffold';
CurveCalcMethod = 'Average';

if FlagSave
    CellColorVec = cell(NumIter,1);
    filename = strcat('C:\Documents\MATLAB\Thesis\3D CCTG\Output data\',ShapeType);
end

Vgrowth = zeros(NumIter,1); % this array stores the filled up volume at each iteration
VGrowthRate = zeros(NumIter,1);
AvertexArray = zeros(NumIter,1); % this array stores the shape surface area at each iteration 
Lambda = zeros(NumIter,1);
TotTime = 0;              % the total time for the simulation

%% Generate the example
[Faces{1}, Vertices{1}, Vtot, NormRatio] = CreateShape(ShapeType);
if FlagDegradation
    temp = Faces{1}(:,2);
    Faces{1}(:,2) = Faces{1}(:,3);
    Faces{1}(:,3) = temp;
end
if FlagPartView
    Faces{1}(1,4) = 0;
    % marking the faces that are included inside a certain area in the mesh
    MinMaxRatio = [0.333 0.666; 0.333 0.666; 0.333 0.666];
%     MinMaxRatio = [0 1; 0.16 0.84; 0 1];
    Faces{1} = DetectInboundMesh(Faces{1}, Vertices{1}, MinMaxRatio);
end
Vertices{1}(1,4:7) = 0;
[emag] = ecalc(Faces{1},Vertices{1}(:,1:3)); 
DeltaRemesh = min(min(emag));
DeltaRemeshMax = 1.3*max(max(emag));
maxis = 1.3*max(max(abs(Vertices{1}(:,1:3))));

%% detect mesh edges for later use
if FlagBoundary
    [Vertices{1}, BoundaryVal, BoundStrVer] = DetectBoundaries(Faces{1},Vertices{1}, 0.5*DeltaRemesh);
end
%% saving the initial geometry for display
if FlagDispInitGeom
    InitGeom = {Faces{1} Vertices{1}};
else
    InitGeom = [];
end
%% open a figure
fig_h = figure('name',ShapeType,'numbertitle','off','color',[0.75 0.75 0.75]);
ax = axes('DataAspectRatio', [1,1,1]);
if FlagVideo
%     set(gcf, 'Position', [0,0,1280,1280]);
    set(gcf, 'Position', [0,0,1280,720]);
    maxis = 1.0*max(max(abs(Vertices{1}(:,1:3))));
else
    set(gcf, 'OuterPosition', get(0, 'Screensize'));
end
axis tight

%% run the algorithm
TotVgrowth = 0;
ii = 1;
while TotTime < TimeLimit
    PlotMesh(Faces{ii}, Vertices{ii}, [], maxis, [], FlagBoundary, FlagPartView, [])
    if FlagBoundary % reflecting the boundaries to calculate the curvatures and normals
        OrigFacesLength = length(Faces{ii}(:,1));
        OrigVertiLength = length(Vertices{ii}(:,1));
        for ll = 1:6
            [RefFaces, RefVertices] = ReflectBoundaryMesh(Faces{ii}, Vertices{ii}, [ceil(0.5*ll) BoundaryVal(ll)], 0.5*DeltaRemesh);
            Faces{ii} = [Faces{ii}; RefFaces];
            Vertices{ii} = [Vertices{ii}; RefVertices];
        end
    end
    PlotMesh(Faces{ii}, Vertices{ii}, [], maxis, [], FlagBoundary, FlagPartView, [])
    
    % calculating the curvatures and normals
    [FaceNormals] = CalcFaceNormals(Faces{ii}(:,1:3),Vertices{ii}(:,1:3));
    [VertexNormals, Avertex, Acorner,up,vp] = CalcVertexNormals(Faces{ii}(:,1:3),Vertices{ii}(:,1:3),FaceNormals);
    [VertexSFM, wfp] = CalcCurvature(Faces{ii}(:,1:3),Vertices{ii}(:,1:3),VertexNormals,FaceNormals,Avertex,Acorner,up,vp);
            
    if FlagBoundary % after calculating the curvatures and normals, cut off the reflected entities
        Faces{ii} = Faces{ii}(1:OrigFacesLength,:);
        Vertices{ii} = Vertices{ii}(1:OrigVertiLength,:);
        VertexNormals = VertexNormals(1:OrigVertiLength,:);
        VertexSFM = VertexSFM(:,:,1:OrigVertiLength);
        Avertex = Avertex(1:OrigVertiLength);
    end
        
    %% Calculating the main curvature in a point and making a step in the direction
    if FlagPartView
        PartVertices = unique(reshape(Faces{ii}(Faces{ii}(:,4)==1,1:3),[numel(Faces{ii}(Faces{ii}(:,4)==1,1:3)),1]));
    else
        PartVertices = [];
    end
    if FlagAnaSol
        [CurvatureVec, ColorVec, AnaSolution] = CurveMethod(CurveCalcMethod, VertexSFM, AnalyticParam, PartVertices);
        if AnaSolution
            break
        end
    else
        [CurvatureVec, ColorVec, AnaSolution] = CurveMethod(CurveCalcMethod, VertexSFM);
    end
%     ColorVec = (squeeze(VertexSFM(1,1,:).*VertexSFM(2,2,:))).';
    if ii == 1
        minmaxColor = [min(ColorVec) max(abs(ColorVec))];
    end
    if FlagDegradation
        Lambda(ii) = LFactor*DeltaRemesh;
        Step = -Lambda(ii)*VertexNormals;
    else
        Lambda(ii) = -LFactor*DeltaRemesh/min(min(CurvatureVec),-max(CurvatureVec)); % setting a value of lambda according to the maximal avg curvature
        Step = Lambda(ii)*VertexNormals.*CurvatureVec;
    end
    Vertices{ii+1} = Vertices{ii} - [Step, zeros(size(Vertices{ii},1),4)];
    
    %% Calculating the total volume filled up in this iteration
    if FlagPartView
        Avertex = CalcVertexAreas(Faces{ii}(Faces{ii}(:,4)==1,1:3),Vertices{ii}(:,1:3));
        AvertexNew = CalcVertexAreas(Faces{ii}(Faces{ii}(:,4)==1,1:3),Vertices{ii+1}(:,1:3));
        Vgrowth(ii) = sum(((Avertex+AvertexNew)/2).*sqrt(sum(Step.^2,2)));
    else
        AvertexNew = CalcVertexAreas(Faces{ii}(:,1:3),Vertices{ii+1}(:,1:3));
        Vgrowth(ii) = sum(((Avertex+AvertexNew)/2).*sqrt(sum(Step.^2,2)));
    end
    AvertexArray(ii) = sum(Avertex);
    TotVgrowth = TotVgrowth + Vgrowth(ii)/(NormRatio^3);
    VGrowthRate(ii) = Vgrowth(ii)/(NormRatio*Lambda(ii));
    if VGrowthRate(ii) < VGrowthRateLim*VGrowthRate(1)
        break
    end
    TotTime = TotTime + Lambda(ii)/(NormRatio^2);
    %% Remeshing the shape to avoid clustering of vertices
    if FlagBoundary
        [Faces{ii+1}, Vertices{ii+1}, BoundStrVer] = Remesh(Faces{ii}, Vertices{ii+1}, DeltaRemesh, DeltaTheta, DeltaRemeshMax, BoundStrVer);
        % imposing special boundary conditions
        [Vertices{ii+1}, BoundStrVer,P] = SpecialBoundaryCond(Vertices{ii+1} ,BoundStrVer);
    else
        [Faces{ii+1}, Vertices{ii+1}] = Remesh(Faces{ii}, Vertices{ii+1}, DeltaRemesh, DeltaTheta, DeltaRemeshMax);
    end  
    %% Plotting the mesh
    if ~mod(ii,1)   
        PlotMesh(Faces{ii}, Vertices{ii}, ColorVec, maxis, minmaxColor, FlagBoundary, FlagPartView, InitGeom)
        % writing video files
        drawnow
        if FlagVideo
            images{ii} = getframe(gcf);
        end
    end
    if FlagSave
        CellColorVec{ii} = ColorVec;       
        file = {Faces, Vertices, CellColorVec, Lambda, Vgrowth};
        save(filename, 'file')
    end
    ii = ii+1;
end
%% Reorganizing the output
Faces = Faces(1:ii-1);
Vertices = Vertices(1:ii-1);
if AnaSolution
    if AnaSolution(1) == 1
        Vvacant = 4*pi*AnaSolution(2)^3/3;
        Lambda(ii) = (AnaSolution(2)^2)/4;
    elseif AnaSolution(1) == 2
        % in the case of a cylinder first the length of the cylinder is
        % estimated using shapes' surface
        if ii ~= 1
            SurfArea = sum(AvertexNew);
        else
            SurfArea = sum(Avertex);
        end
        Vvacant = SurfArea*AnaSolution(2)/2;
        Lambda(ii) = (AnaSolution(2)^2)/2;
    end
    Vgrowth(ii) = Vvacant;
    Vgrowth = Vgrowth(1:ii);
    AvertexArray = AvertexArray(1:ii);
    Lambda = Lambda(1:ii);
else
    Vgrowth = Vgrowth(1:ii-1);
    AvertexArray = AvertexArray(1:ii-1);
    Lambda = Lambda(1:ii-1);
end

%% producing a video
if FlagVideo
    images = images(1:ii-1);
    writerObj = ProduceVideo(images,ShapeType,Lambda,TotVideoLength);
end

%% calculating the loss function
if FlagDegradation
    TotalLoss = Vtot/TotVgrowth-1;
else
    Vgrowth = Vgrowth/(NormRatio^3);
    AvertexArray = AvertexArray/(NormRatio^2);
    Lambda = Lambda/(TimeFactor*NormRatio^2);
    TrueTotTime = sum(Lambda);
    TotalLoss = TrueTotTime*(Vtot/TotVgrowth);
    disp('Total true time to fill:'); 
    disp(num2str(TrueTotTime));
end
%% saving the data
if FlagSave   
    file = {Faces, Vertices, CellColorVec, Lambda, Vgrowth, Vtot, NormRatio};
    save(filename, 'file')
end

% Surprise_Motherfucker;
