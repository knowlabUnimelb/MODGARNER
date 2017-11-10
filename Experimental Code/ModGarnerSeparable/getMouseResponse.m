function [charx, t] = getMouseResponse(screenparms)

FlushEvents('mouseDown');
startSecs = GetSecs;
while 1
    [x, y, buttons] = GetMouse(screenparms.window);
    %     if any(buttons(1:2))
    %         charx = find(buttons(1:2));
    %         t = GetSecs - startSecs;
    %         break;
    %     end
    if any(buttons([1 3]))
        charx = find(buttons([1 3]));
        t = GetSecs - startSecs;
        break;
    end
end

FlushEvents('keyDown');