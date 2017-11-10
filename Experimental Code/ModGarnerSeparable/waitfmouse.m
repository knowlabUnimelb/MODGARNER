function waitfmouse(screenparms)
global ptb3


while KbCheck; end
FlushEvents('mouseDown')
showtext(screenparms, 25, 'Press a button on the mouse to continue', 1, round(screenparms.rect(4)/2.5), -round(screenparms.rect(3)/3)); if ptb3; Screen('Flip', screenparms.window); end

buttons = [1 1 1];
while 1
    while any(buttons)
        [x,y,buttons] = GetMouse(screenparms.window);
    end
    [x,y,buttons] = GetMouse(screenparms.window);
    if sum(double(buttons)) > 0
        break
    end
end
if ptb3; Screen('FillRect', screenparms.window, screenparms.screen); else Screen('FillRect'); end
FlushEvents('mouseDown')
while KbCheck; end