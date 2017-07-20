function [ padMatrix ] = Mpad( oMatrix, sizeEX, sizeEY )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
oMatrix = sparse(oMatrix);

[sizeX, sizeY] = size(oMatrix);

zVMatrix = sparse(sizeEY, sizeX);
zHMatrix = sparse(sizeY, sizeEX);
zXMatrix = sparse(sizeEX, sizeEY);

padMatrix = [oMatrix; zVMatrix];

z2Matrix = [zHMatrix;zXMatrix];

padMatrix = [padMatrix z2Matrix];

end

