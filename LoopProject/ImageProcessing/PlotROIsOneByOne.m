% Plot one ROI and wait for keyboard inputs to plots the next.
% Keyboard input one must be y/n and 2 must be a number btwn
% 0-9, to indicate the level of confidence in identifying the beads


% for i = 1:data.Info.nROIs
%     if data.s2p.iscell(i)
%         figure;hold on
%         for ii = 1 : size(data.s2p.stat{1, 1}.xpix,2)
%             imagesc(double(min(data.s2p.stat{1, i}.xpix(1,ii))),...
%                 double(min(data.s2p.stat{1, i}.ypix(1,ii))),...
%                 data.s2p.meanImg(double(min(data.s2p.stat{1, i}.xpix(1,ii))),...
%                 double(min(data.s2p.stat{1, i}.ypix(1,ii)))))
%         end
%     end
% end

% beads = NaN(sum(data.s2p.iscell),3);
counter = 1;
for i = 4:data.Info.nROIs
    if data.s2p.iscell(i)
        beads(counter,1) = counter;
        h1 = figure;
        subtightplot(1,2,1);hold on
        imagesc(data.s2p.meanImg)
        x = double(data.s2p.stat{1, i}.xpix);
        y = double(data.s2p.stat{1, i}.ypix);
                plot(x,y,'r');
        axis equal
        xlim([0 size(data.s2p.meanImg,1)]); ylim([0 size(data.s2p.meanImg,2)]);
        axis ij
        colormap('gray')
        subtightplot(1,2,2);hold on
        imagesc(data.s2p.meanImg)
        text(data.s2p.stat{1, i}.med(2),data.s2p.stat{1, i}.med(1),num2str(i),...
            'Color',[1 0 0],'HorizontalAlignment','center','FontSize',12);
        axis equal
        xlim([0 size(data.s2p.meanImg,1)]); ylim([0 size(data.s2p.meanImg,2)]);
        axis ij
        colormap('gray')
        set(gcf,'Units','Normalized','Position',[0.35 0.32 0.65 0.6])
        k = waitforbuttonpress;
        % 121 y
        % 110 n
        value = double(get(gcf,'CurrentCharacter'));
        if value == 121
            beads(counter,2) = true;
        elseif value == 110
            beads(counter,2) = false;
        end
        k = waitforbuttonpress;
        confidence = double(get(gcf,'CurrentCharacter'));
        beads(counter,3) = (confidence-48)/10;
        fprintf('beads: %g, confidence: %g \n', beads(counter,2), beads(counter,3))
        close(h1);
        counter = counter+1;
    end
end