function [Faces] = DetectInboundMesh(Faces, Vertices, MinMaxRatio)
% This function marks the faces that are located inside a
% certain area.
% MinMaxVal = [Minimal x Value, Maximal x Value
%              Minimal y Value, Maximal y Value
%              Minimal z Value, Maximal z Value]
eps = 2*10^-2;
MinMaxVal = [min(Vertices(:,1))-eps max(Vertices(:,1))+eps
             min(Vertices(:,2))-eps max(Vertices(:,2))+eps
             min(Vertices(:,3))-eps max(Vertices(:,3))+eps];
Lengths = MinMaxVal(:,2) - MinMaxVal(:,1);
MinMaxVal = MinMaxVal(:,1) + Lengths.*MinMaxRatio;

IncVerticesNum = zeros(length(Vertices(:,1)),1);
jj=1;
for ii=1:length(Vertices(:,1))
    if Vertices(ii,1)>=MinMaxVal(1,1)&&Vertices(ii,1)<=MinMaxVal(1,2)&&Vertices(ii,2)>=MinMaxVal(2,1)&&Vertices(ii,2)<=MinMaxVal(2,2)&&Vertices(ii,3)>=MinMaxVal(3,1)&&Vertices(ii,3)<=MinMaxVal(3,2)
        IncVerticesNum(jj) = ii;
        jj=jj+1;
    end
end
IncVerticesNum = IncVerticesNum(1:jj-1);

for ii=1:length(Faces(:,1))
    if sum(Faces(ii,1)==IncVerticesNum) + sum(Faces(ii,2)==IncVerticesNum) + sum(Faces(ii,3)==IncVerticesNum)==3
        Faces(ii,4) = 1;
    end
end

