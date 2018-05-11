function [data, parms] = generateTestData(model, nd, stimloc, varargin)

optargs = {1.5, 1.5, 1, 1, .75, .25, .25, .1, .3, .5, 1, .5, .75, .75};
newVals = cellfun(@(x) ~isempty(x), varargin); % skip any new inputs if they are empty
optargs(newVals) = varargin(newVals); % now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
[db1, db2, sp1, sp2, A, bMa1, bMa2, s, t0, pX, m, pSer, Aser, Apar] = optargs{:}; % Place optional args in memorable variable names

parms = struct('db1', db1, 'db2', db2, 'sp1', sp1, 'sp2', sp2,...
               'A', A, 'bMa1', bMa1, 'bMa2', bMa2, 's', s, 't0', t0,...
               'pX', pX, 'm', m, 'pSer', pSer, 'Aser', Aser, 'Apar', Apar);

%% Parameters
stoppingrule = model(end-1:end);
db = [parms.db1, parms.db2];
sp = [parms.sp1, parms.sp2];

flag = false;
if strcmp(model, 'coactive  ')
    flag = true;
end
vc = grt_dbt(stimloc, db, sp, flag); % Drift rate for the "A" accumulator

if strcmp(model, 'mixedSPst')
    servc = vc;
    parvc = grt_dbt(stimloc, db, parms.m*sp, flag); % Drift rate for the "A" accumulator 
end

%% Simulate data
n.items = size(vc,1);

data = struct('rt', [], 'resp', [], 'item', [], 'stimloc', stimloc);
for i = 1:n.items
    if strcmp(model, 'serialst')
        simparms = struct('vc', vc(i,:), 'A', parms.A,...
            'bMa', [parms.bMa1, parms.bMa2], 's', parms.s, 't0', parms.t0, 'pX', parms.pX);
    elseif strcmp(model, 'mixedSPst')
        simparms = struct('servc', servc(i,:), 'parvc', parvc(i,:), 'Aser', parms.Aser, 'Apar', parms.Apar,...
            'bMa', [parms.bMa1, parms.bMa2], 's', parms.s, 't0', parms.t0, 'pX', parms.pX, 'pSer', parms.pSer);
    else
        simparms = struct('vc', vc(i,:), 'A', parms.A, 'bMa', [parms.bMa1, parms.bMa2], 's', parms.s, 't0', parms.t0);
    end

    eval(sprintf('[rt, resp] = sim%s(nd, simparms, stoppingrule);', model(1:end-2)))
    
    data.rt = [data.rt; rt];
    data.resp = [data.resp; resp];
    data.item = [data.item; i * ones(nd, 1)];
end

maxrt = 10;
data.resp(data.rt > maxrt) = [];
data.item(data.rt > maxrt) = [];
data.rt(data.rt   > maxrt) = [];
