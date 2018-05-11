function summarizePosteriorsFormatted(phi, n)

names = {'c', 'cf', 'bh', 'A', 'bMa1', 'bMa2', 'vs', 't0',...
    'wp', 'b1', 'b2', 'b3', 'b4', 'b5'};

table = cell(numel(names), 4);
cnt = 1;
for i = 1:numel(names)

    switch names{cnt}
        case {'c', 'cf', 'bh', 'A', 'bMa1', 'bMa2', 'vs', 't0'}
            temp1 = phi.(sprintf('%s_mu',names{cnt}))(:,n.burnin:end);
            temp1 = temp1(:);
            
            temp2 = phi.(sprintf('%s_sigma',names{cnt}))(:,n.burnin:end);
            temp2 = temp2(:);

            
            fprintf('m = [%3.2f, %3.2f], s = [%3.2f, %3.2f]\n',...
                round(prctile(temp1, [2.5]), 2), round(prctile(temp1, [97.5]), 2),...
                round(prctile(temp2, [2.5]), 2), round(prctile(temp2, [97.5]), 2))
        case {'w', 'wc', 'wp', 'b1', 'b2', 'b3', 'b4', 'b5'}
            temp1 = phi.(sprintf('%s_a',names{cnt}))(:,n.burnin:end);
            temp1 = temp1(:);
            
            temp2 = phi.(sprintf('%s_b',names{cnt}))(:,n.burnin:end);
            temp2 = temp2(:);
            
            fprintf('a = [%3.2f, %3.2f], b = [%3.2f, %3.2f]\n',...
                round(prctile(temp1, [2.5]), 2), round(prctile(temp1, [97.5]), 2),...
                round(prctile(temp2, [2.5]), 2), round(prctile(temp2, [97.5]), 2))            
    end
    cnt = cnt + 1; 
    
end
