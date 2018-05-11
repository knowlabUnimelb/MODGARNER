function flag = checkParms(parms, check, names)

if any(check < 0)
    flag = true;
elseif any(strcmp(names, 'w')) && any(parms.('w') > 1)
    flag = true;
elseif any(strcmp(names, 'wc')) && any(parms.('wc') > 1)
    flag = true;
elseif any(strcmp(names, 'b1')) && any(parms.('b1') > 1)
    flag = true;
elseif any(strcmp(names, 'b2')) && any(parms.('b2') > 1)
    flag = true;
elseif any(strcmp(names, 'b3')) && any(parms.('b3') > 1)
    flag = true;
elseif any(strcmp(names, 'b4')) && any(parms.('b4') > 1)
    flag = true;
elseif any(strcmp(names, 'b5')) && any(parms.('b5') > 1)
    flag = true;
elseif any(strcmp(names, 'wp')) && any(parms.('wp') > 1)
    flag = true;
else
    flag = false;
end