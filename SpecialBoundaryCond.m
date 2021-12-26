function [Vertices, BoundStrVer, P] = SpecialBoundaryCond(Vertices, BoundStrVer)
% this function checks for special boundary conditions on the function and
% adjusts boundaries/locations if necessary.
%% checking for vanishing boundaries
ii=1;
while ii <= length(BoundStrVer)
    if length(BoundStrVer{ii}) < 4
        Vertices(BoundStrVer{ii},4:6) = 0;
        BoundStrVer = [BoundStrVer(1:ii-1); BoundStrVer(ii+1:end)];
        continue
    end
    ii=ii+1;
end
%% manually setting the boundary vertices which are also corners
CornerFlag = 0;
if CornerFlag
    CornVerInd = find(sum(Vertices(:,4:6),2)==2);
    P = zeros(length(CornVerInd),3);
    for ii = 1:length(CornVerInd)
        Neighbors = Find4Nieghbors(CornVerInd(ii), BoundStrVer);   
        % second order estimation using the four closest neighbors
        % this according to the formula: O = A+AB*AP/|AB,AB|
        % P is an estimation based on the four neigbors
        Dist= [sqrt(sum((Vertices(Neighbors(2),1:3)-Vertices(Neighbors(1),1:3)).^2)) sqrt(sum((Vertices(Neighbors(4),1:3)-Vertices(Neighbors(5),1:3)).^2))];
        AvgDistance = sum(Dist(:))/2;
        Points = [Vertices(Neighbors(2),1:3) + AvgDistance*(Vertices(Neighbors(2),1:3)-Vertices(Neighbors(1),1:3))/Dist(1)
                  Vertices(Neighbors(4),1:3) + AvgDistance*(Vertices(Neighbors(4),1:3)-Vertices(Neighbors(5),1:3))/Dist(2)];
        P(ii,:) = sum(Points)/2;
        % first order estimation using the two closest neighbors
        % using the formula: O = A+AB*AP/|AB,AB|
        % P is the average location of the two neighbors
%         P = (Vertices(Neighbors(2),1:3)+Vertices(Neighbors(4),1:3))/2;
    % projecting the point upon the corner edge
%         A = Vertices(Neighbors(3),1:3);
%         AB = A.*(~Vertices(Neighbors(3),4:6));
%         normAB = AB/dot(AB,AB);
%         Vertices(CornVerInd(ii),1:3) = A+dot(AB,P(ii,:)-A)*normAB; 
    end
else
    P = [];
end
end

function [Neighbors] = Find4Nieghbors(VerInd, BoundStrVer)
for ii = 1:length(BoundStrVer)
    Index = find(BoundStrVer{ii}(2:end) == VerInd);
    if Index
        Neighbors = BoundStrVer{ii}(mod(Index-2:Index+2,length(BoundStrVer{ii})-1)+1);
    end
end
end
