function showInstructions(screenparms, instructionFile, advance, waittime)
if nargin == 2
    advance = 'space';
end

imageInfo = imfinfo(instructionFile);
instructionImage = imread(instructionFile, imageInfo.Format);
if size(instructionImage, 3) == 1
    instructionImage = ind2rgb(imread(instructionFile, imageInfo.Format), imageInfo.Colormap) * 255;
end

if imageInfo.Height<0
    imageInfo.Height=-imageInfo.Height;
end

if imageInfo.Width<0
    imageInfo.Width=-imageInfo.Width;
end

[imageRect,dh,dv] = CenterRect([0 0 imageInfo.Width imageInfo.Height], screenparms.rect);
if screenparms.rect(3)/imageInfo.Width < 1 || screenparms.rect(4)/imageInfo.Height < 1
    newWidth = imageInfo.Width * min([screenparms.rect(3)/imageInfo.Width, screenparms.rect(4)/imageInfo.Height]);
    newHeight = imageInfo.Height * min([screenparms.rect(3)/imageInfo.Width, screenparms.rect(4)/imageInfo.Height]);
    [imageRect,dh,dv] = CenterRect([0 0 newWidth newHeight], screenparms.rect);
end



instructionTexture = Screen('MakeTexture', screenparms.window, instructionImage);
Screen('DrawTexture', screenparms.window, instructionTexture, [], imageRect);
Screen('Flip', screenparms.window);

switch advance
    case 'space'
        try
            pressSpace; Screen('Flip', screenparms.window); 
        catch
            pause; FillScreen(screenparms); Screen('Flip', screenparms.window); 
        end
    case 'time'
        if nargin == 3
            waittime = 1;
        end
        WaitSecs(waittime);
    case 'RTBox'
          RTBox('clear'); % clear buffer and sync clocks before stimulus onset
          while ~any(RTBox('ButtonDown')); WaitSecs(0.01); end; 
          WaitSecs(.1);
    case 'PsychRTBox'
%           PsychRTBox('clear', rthandle); % clear buffer and sync clocks before stimulus onset
          while ~any(PsychRTBox('ButtonDown')); WaitSecs(0.01); end; 
          WaitSecs(.1);          
end
Screen('Close', instructionTexture); 