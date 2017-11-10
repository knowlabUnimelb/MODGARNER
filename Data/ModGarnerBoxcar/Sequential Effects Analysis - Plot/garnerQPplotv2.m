function garnerQPplotv2(mq, task, item, str)
% v2 works with plotSeqData..._v2

if task < 3
    iq = mq{task}(mq{task}(:,2) == item, :);
    for i = 1:5
        plot(1:6, iq(:,3+i-1), str)
        hold on
    end
    set(gca,'XTickLabel', num2str(iq(:,1)), 'XLim', [.5 6.5])
    xlabel('Previous Item', 'FontSize', 14)
    ylabel('RT (msec)', 'FontSize', 14)
else
    sets = [1:2:12; 2:2:12];

    linesettings = str; %{'-ok','-^r'};
    offsets = [-.15, .15];
    for ss = 1:size(sets,1)
        iq = mq{task}(mq{task}(:,2) == item & ismember(mq{task}(:,1), sets(ss,:)), :);
        for i = 1:5
            if ss == 1
                plot((1:6)+offsets(ss), iq(:,3+i-1), linesettings{ss})
            else
                plot((1:6)+offsets(ss), iq(:,3+i-1), linesettings{ss}, 'MarkerFaceColor', linesettings{ss}(end))
            end
            hold on
        end
        set(gca,'XTickLabel', num2str(iq(:,1)), 'XLim', [.5 6.5])
        xlabel('Previous Item', 'FontSize', 14)
        ylabel('RT (msec)', 'FontSize', 14)
    end
end