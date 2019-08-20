function [ind,e] = makeWin(data, range)

ind = data > range(1) & data < range(2); %Returns the indices of the specified range
e = data(ind); %Returns the values within that range