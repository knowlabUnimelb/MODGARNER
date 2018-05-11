function idx = findMulti(x, y)

for i = 1:numel(y)
   idx(i) = find(x <= y(i), 1, 'last'); 
end