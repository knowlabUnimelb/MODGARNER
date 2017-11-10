function waitfkey (screenparms, stimsize, msg)
global ptb3

while KbCheck; end
temp = 100;%expinfo.diameter; % temporarily alter font size
stimsize = 100/5;%expinfo.diameter/5;
showtext(screenparms, stimsize, msg, 1, 0,-100); if ptb3; Screen('Flip', screenparms.window); end
stimsize = temp;
while 1
	[keyIsDown,timeSecs,keyCode] = KbCheck;
	if keyIsDown
		c = KbName(keyCode);
        while KbCheck; end
        if strcmp(c,'space')
            break;
        end
	end
end
if ptb3; Screen('FillRect', screenparms.window, screenparms.screen); else Screen('FillRect'); end
while KbCheck; end