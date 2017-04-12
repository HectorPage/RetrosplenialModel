function plotRates(seconds,tsize,rates, CellType, NumCells, save_mode, test_type)
%plotRates - Plots rates from various regions
% Created by Dr. Hector JI Page 24/10/16
%modified by Dr. Hector JI Page 03/02/2017

%Test type tells the model which sort of test has been done on the data to
%plot:
%1. SO_test = no ADN PI input
%2. disorientation = initial ADN packet in random position
%3. none = no test done

test_type = test_type(1); %only need first letter

x = tsize:(tsize*100):seconds; %Used for models saving every 100th timestep
%xtick = [0,steps/4, steps/2, (steps/4)*3, steps];
xtick = [0,seconds/4, seconds/2, (seconds/4)*3, seconds];

y = 1:NumCells;
[X, Y] = meshgrid(x, y);

figure();
h1 = surf(X, Y, rates);

set(h1, 'LineStyle', 'none');
xlabel('Time (s)', 'Fontsize', 32);
set(gca, 'XTick', xtick, 'Fontsize', 24);
%set(gca, 'XTickLabel',xticklabel,'Fontsize',24);
ylabel('Cell', 'Fontsize', 32);
ylim([1 NumCells]);
xlim([0 seconds]);
colormap(1-gray);
view(0, 90);
title([CellType,' Cell Rates'], 'FontSize', 32);

if (strcmpi(save_mode, 'save'))
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    if strcmpi(test_type,'s')
        filename = ['_',CellType,'Rates_Test'];
    elseif strcmpi(test_type,'d')
        filename = ['_',CellType,'Rates_disor'];
    else
        filename = ['_',CellType,'Rates'];
    end
    saveas(gcf,filename, 'fig');
    close(gcf);
end



end

