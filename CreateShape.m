function [Faces, Vertices, Vtot, NormRatio] = CreateShape(Shape)
% this function creates a mesh of a certain shape and returns it
% Vvacant - the empty volume (sometimes taken from a CAD model)
% SurfArea - the surface area of the shape
switch Shape
    case 'Peaks'
        %% Generate a peaks surface
        StepSize=0.5;
        GridSize=20;
        Grid = 1:StepSize:GridSize;
        [x,y] = meshgrid(Grid,Grid);
        Faces = delaunay(x,y);
        z = peaks(size(Grid,2));
        Vertices =[x(1:end)',y(1:end)',z(1:end)'];
        Vtot = 1;
        NormRatio = 1;
    case {'Sphere', 'Ellipsoid', 'Tetrahedron', 'Cube', 'Octahedron', 'Dodecahedron', 'Icosahedron', 'Stellated Dodecahedron', 'Stellated Octahedron'}
        %% One of the 3D model surfaces/polyhedrons
        [Faces, Vertices, ~] = Read_Data(Shape);
        temp = Faces(:,2);
        Faces(:,2) = Faces(:,3);
        Faces(:,3) = temp;
        
        Vtot = 1;
        SurfArea = 4*10^-3;
        [Vertices, NormRatio] = NormPoints(Vertices);
    case 'Regular_Scaffold'
        [Faces, Vertices, ~] = Read_Data(Shape);
        
        Vtot = (1.6*10^-4)*(2*10^-4)^2;
        SurfArea = 1;
        [Vertices, NormRatio] = NormPoints(Vertices);
%         TotVideoLength = 14.8;
%         TimeLimit = 4e-09;;
%         Lfac = 0.14;
    case {'Lattice_Scaffold_Degradation', 'Lattice_Scaffold_Degradation_2'}
        [Faces, Vertices, ~] = Read_Data(Shape);
        
        Vtot = 3028.89*10^-9;   % the initial scaffold volume [m^3]
        SurfArea = 1;
        [Vertices, NormRatio] = NormPoints(Vertices);
%         TotVideoLength = ?;
%         TimeLimit = 4.5e-05;
%         LFactor = 0.3;
    case {'Square_Scaffold_Degradation', 'Square_Scaffold_Degradation_2'}
        [Faces, Vertices, ~] = Read_Data(Shape);
        
        Vtot = 1.4336e-10;   % the initial scaffold volume [m^3]
        SurfArea = 1;
        [Vertices, NormRatio] = NormPoints(Vertices);
%         TotVideoLength = ?;
%         TimeLimit = 4.1992e-08;
%         LFactor = 0.3;
    case 'BCC_Scaffold_Degradation'
        [Faces, Vertices, ~] = Read_Data(Shape);
        
        Vtot = 1.9919e-7;   % the initial scaffold volume [m^3]
        SurfArea = 1;
        [Vertices, NormRatio] = NormPoints(Vertices);
%         TotVideoLength = ?;
%         TimeLimit = 4.5e-05;
%         LFactor = 0.3;
    case {'Cylindrical_Scaffold_Degradation', 'Cylindrical_Scaffold_Degradation_2'}
        [Faces, Vertices, ~] = Read_Data(Shape);
        
        Vtot = 4.5724e-7;   % the initial scaffold volume [m^3]
        SurfArea = 1;
        [Vertices, NormRatio] = NormPoints(Vertices);
%         TotVideoLength = ?;
%         TimeLimit = 9.3646e-6;
%         LFactor = 0.3;
    case 'Scaffold_Unit_Plus_Plus'
        Faces = xlsread('Scaffold_Unit_Plus_Plus_Elements_Data','C:E');
        Vertices = xlsread('Scaffold_Unit_Plus_Plus_Nodes_Data','B:D');
        
        Vtot = (1.6*10^-4)*(1.9635*10^-4)^2;
        SurfArea = 9.957*10^-8;
        [Vertices, NormRatio] = NormPoints(Vertices);
%         TotVideoLength = 14.8;
%         TimeLimit = 0.39;
%         Lfac = 0.14;
%         MinMaxRatio = [0 1; 0.16 0.84; 0 1];
    case 'Scaffold_4_Units'
        Faces = xlsread('Scaffold_4_Units_Elements_Data','C:E');
        Vertices = xlsread('Scaffold_4_Units_Nodes_Data','B:D');
        
        Vtot = (1.6*10^-4)*(1.9635*10^-4)^2;
        SurfArea = 1.5485*10^-7;
        [Vertices, NormRatio] = NormPoints(Vertices);
case 'Scaffold_27_Units'
        Faces = xlsread('Scaffold_27_Units_Elements_Data','C:E');
        Vertices = xlsread('Scaffold_27_Units_Nodes_Data','B:D');
        
        Vtot = (1.6*10^-4)*(1.9635*10^-4)^2;
        SurfArea = 9.957*10^-8;
        [Vertices, NormRatio] = NormPoints(Vertices);
%         TimeLimit = ;
    case 'BCC_Scaffold'
        [Faces, Vertices, ~] = Read_Data(Shape);
        
        Vtot = 0.01^3;
        SurfArea = 1;
        [Vertices, NormRatio] = NormPoints(Vertices);
%         TotVideoLength = ;
%         TimeLimit = 1.85*10^-5;
%         Lfac = ;
    case 'Duck'
        %% Generate a Box surface
        Faces = xlsread('Duck_Elements_Data','C:E');
        Vertices = xlsread('Duck_Nodes_Data','B:D');
        SurfArea = 1;
        Vtot = 1;
        [Vertices, NormRatio] = NormPoints(Vertices);
    case 'Bannana'
        %% Generate a Box surface
        Faces = xlsread('Bannana_Elements_Data','C:E');
        % adjusting the second and third column such that normal will point inside
        temp = Faces(:,2);
        Faces(:,2) = Faces(:,3);
        Faces(:,3) = temp;
        Vertices = xlsread('Bannana_Nodes_Data','B:D');
        
        Vtot = 1;
        SurfArea = 1;
        [Vertices, NormRatio] = NormPoints(Vertices);
    case 'Vase'
        %% Generate a Vase surface
        Faces = xlsread('Vase_Elements_Data','C:E');
        % adjusting the second and third column such that normal will point inside
        temp = Faces(:,2);
        Faces(:,2) = Faces(:,3);
        Faces(:,3) = temp;
        
        Vertices = xlsread('Vase_Nodes_Data','B:D');
        Vertices = Vertices - mean(Vertices);
        
        Vtot = 1;
        SurfArea = 1;
        [Vertices, NormRatio] = NormPoints(Vertices);
%         RecTime = 0.06;
%         TotVideoLength = 18.1;
%         view([0 0.75 -1])
    case 'Infinity'
        %% Generate a The infinity symbol surface
        Faces = xlsread('Infinity_Elements_Data','C:E');
        % adjusting the second and third column such that normal will point inside
        temp = Faces(:,2);
        Faces(:,2) = Faces(:,3);
        Faces(:,3) = temp;
        
        Vertices = xlsread('Infinity_Nodes_Data','B:D');
        SurfArea = 1;
        [Vertices, NormRatio] = NormPoints(Vertices);
%         RecTime = 0.009;
%         TotVideoLength = 9;
%         caxis([-0.38 0.49])
    case 'Triangle'
        %% Generate a 2D test of a triangle
        Faces = xlsread('Triangle_Elements_Data','C:E');
        % adjusting the second and third column such that normal will point inside
        temp = Faces(:,2);
        Faces(:,2) = Faces(:,3);
        Faces(:,3) = temp;
        Vertices = xlsread('Triangle_Nodes_Data','B:D');
        SurfArea = 4*10^-3;
        Vtot = 1;
        [Vertices, NormRatio] = NormPoints(Vertices);
    case 'Square'
        %% Generate a 2D test of a square
        Faces = xlsread('Square_Elements_Data','C:E');
        Vertices = xlsread('Square_Nodes_Data','B:D');
        % adjusting the second and third column such that normal will point inside
        temp = Faces(:,2);
        Faces(:,2) = Faces(:,3);
        Faces(:,3) = temp;
        
        Vtot = 1;
        SurfArea = 4*10^-3;
        [Vertices, NormRatio] = NormPoints(Vertices);
    case 'Hexagon'
        %% Generate a 2D test of a hexagon
        Faces = xlsread('Hexagon_Elements_Data','C:E');
        Vertices = xlsread('Hexagon_Nodes_Data','B:D');
        % adjusting the second and third column such that normal will point inside
        temp = Faces(:,2);
        Faces(:,2) = Faces(:,3);
        Faces(:,3) = temp;
        
        Vtot = 1;
        SurfArea = 4*10^-3;
        [Vertices, NormRatio] = NormPoints(Vertices);
    case 'Circle'
        %% Generate a 2D test of a Circle
        Faces = xlsread('Circle_Elements_Data','C:E');
        Vertices = xlsread('Circle_Nodes_Data','B:D');
        % adjusting the second and third column such that normal will point inside
        temp = Faces(:,2);
        Faces(:,2) = Faces(:,3);
        Faces(:,3) = temp;
        
        Vtot = 12732.39*10^-9;
        SurfArea = 4*10^-3;
        [Vertices, NormRatio] = NormPoints(Vertices);
    case 'Cross'
        %% Generate a 2D test of a Cross
        [Faces, Vertices, ~] = Read_Data(Shape);
        temp = Faces(:,2);
        Faces(:,2) = Faces(:,3);
        Faces(:,3) = temp;
        
        Vtot = 1;
        SurfArea = 4*10^-3;
        [Vertices, NormRatio] = NormPoints(Vertices);
    case 'Star'
        %% Generate a 2D test of a Star
        [Faces, Vertices, Flag] = Read_Data(Shape);
        temp = Faces(:,2);
        Faces(:,2) = Faces(:,3);
        Faces(:,3) = temp;
        
        Vtot = 1;
        SurfArea = 4*10^-3;
        [Vertices, NormRatio] = NormPoints(Vertices);
     case 'LSM_Scaffold'
        %% Generate a 2D test of a Star
        [Faces, Vertices, Flag] = Read_Data(Shape);
%         temp = Faces(:,2);
%         Faces(:,2) = Faces(:,3);
%         Faces(:,3) = temp;
        
%         TimeLimit = 0.0022        
        Vtot = 0.1^3;
        SurfArea = 4*10^-3;
        [Vertices, NormRatio] = NormPoints(Vertices);
    case 'Bartolo_Comparison'
        Faces = xlsread('Bartolo_Elements_Data','C:E');
        Vertices = xlsread('Bartolo_Nodes_Data','B:D');
        Vertices = Vertices / 10^2;
        % adjusting the second and third column such that normal will point inside
%         temp = Faces(:,2);
%         Faces(:,2) = Faces(:,3);
%         Faces(:,3) = temp;
        
        Vtot = 1;
        SurfArea = 1;
        [Vertices, NormRatio] = NormPoints(Vertices);
    case 'Small_square_scaffold'
        Faces = xlsread('Small_Square_Scaffold_Elements_2','C:E');
        Vertices = xlsread('Small_Square_Scaffold_Nodes_2','B:D');
        
        Vtot = 1;
        SurfArea = 1;
        [Vertices, NormRatio] = NormPoints(Vertices);
    case 'Big_square_scaffold'
        Faces = xlsread('Big_Square_Scaffold_Elements_2','C:E');
        Vertices = xlsread('Big_Square_Scaffold_Nodes_2','B:D');
        
        Vtot = 1;
        SurfArea = 1;
        [Vertices, NormRatio] = NormPoints(Vertices);
    otherwise
        [Faces, Vertices, ~] = Read_Data(Shape);
        [Vertices, NormRatio] = NormPoints(Vertices);
        Vtot = 1;
end
end

function [Vertices, NormRatio] = NormPoints(Vertices)
% this function normalizes a vetex array to a unit sphere
% NormRatio - the normalization ratio (
Vertices = Vertices - mean(Vertices); % translating the array so surround the origin
NormRatio = 1/(max(max(abs(Vertices))));
Vertices = Vertices.*NormRatio;
end


function [Faces, Vertices, Flag] = Read_Data(Shape)
% this function imports a .dat file of a certain mesh and returns:
% Faces - the element array of vertices triangles
% Vertices - the array of node coordinates
% NormRatio - a normalizing ratio used to speed calculations
% Flag - a signal telling whether the input mesh is of acceptable form
% (Flag = 0 usually means a corrupt geometry)
ShapeName = sprintf('%s.dat',Shape); 
Flag = 1;
InitData = readtable(ShapeName,'Delimiter','\n','ReadVariableNames',false); 
idx = ismember(InitData.Var1, '-1');
[idxnum,~] = find(idx);

opts = delimitedTextImportOptions(...
        'Delimiter', ' ',...
        'ConsecutiveDelimitersRule', 'join','LeadingDelimitersRule','ignore',...
        'MissingRule', 'omitvar','LineEnding','\n','ImportErrorRule','omitvar');

InitData = readtable(ShapeName,opts);
% in the cases where there are multipile element tables (idxnum>2)
% or volumetric elements (width(InitData)>16) send failure signal.
if length(idxnum) > 2 || width(InitData) > 16
    Flag = 0;
end
% extracting the Elements and Node tables
NoTable = InitData(26:idxnum(1)-1,2:4);
ElTable = InitData(idxnum(1)+7:idxnum(2)-1,12:15);
%% converting the tables to arrays
Nodes = table2array(NoTable);
Vertices = str2double(Nodes);
Elements = table2array(ElTable);
Faces = str2double(Elements);
% in case of quadrilateral elements send failure signal
if ~isempty(find(Faces(:,3)~=Faces(:,4),1))
    Flag = 0;
end
Faces = Faces(:,1:3);
end