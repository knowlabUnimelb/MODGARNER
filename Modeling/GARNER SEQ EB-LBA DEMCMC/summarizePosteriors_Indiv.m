function table = summarizePosteriors(theta, n, names)

getPostBurn = @(x)(x(:,n.burnin:end));



% table = cell(numel(names), 4);
table = [];
for i = 1:n.subLevParms
    
%     temp = theta.(names{i})(:,n.burnin:end);
%     temp = temp(:);
    
%     table{i,1} = names{i};
%     table{i,2:n.subjects+1} =
    table = [table; cellfun(@(x)(mode(x(:))), cellfun(@(x)getPostBurn(x), theta.(names{i}), 'UniformOutput', false))];
%     table{i,3} = prctile(temp, [2.5]);
%     table{i,4} = prctile(temp, [97.5]);
end
% disp(table)