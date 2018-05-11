function summarizePosteriors(phi, n, names)

table = cell(numel(names), 4);
for i = 1:n.hyperParms
    
    temp = phi.(names{i})(:,n.burnin:end);
    temp = temp(:);
    
    table{i,1} = names{i};
%     table{i,2} = mode(temp);
    table{i,3} = round(prctile(temp, [2.5]), 2);
    table{i,4} = round(prctile(temp, [97.5]), 2);
%     table{i,5} = median(temp);
%     table{i,6} = mean(temp);
end
disp(table)