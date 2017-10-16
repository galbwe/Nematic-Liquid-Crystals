%given 2d array of 1's and 0's, finds the coordinates of all
%entries that are 1
function [X,Y] = findOnes(a)
    X = [];Y = [];
    [N,M] = size(a);
    for row = 1:N
       cols = find(a(row,:) == 1);
       X = [X,cols];
       Y = [Y,row.*ones(1,length(cols))];
    end
end