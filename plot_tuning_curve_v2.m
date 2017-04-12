function plot_tuning_curve_v2(angles,compartments,odour_context,save_mode,cell_number,test_type, celltype,varargin)
%Plots tuning curve generated from model data in the RSC multi-box
%simulations


%Created by and copyright held by Dr. Hector JI Page 24/10/16
%modified by Dr. Hector JI Page 03/02/2017

test_type = test_type(1);

%Test type tells the model which sort of test has been done on the data to
%plot:
%1. SO_test = no ADN PI input
%2. disorientation = initial ADN packet in random position
%3. none = no test done


f = figure();
polar(0,1,'-k'); %setting radius to 1
hold on;
if(compartments>1.5)
    h1 = polar(angles',varargin{1});
    h2 = polar(angles',varargin{2});
else
    h1 = polar(angles',varargin{1});
end

%setting axis linewidth
hlines = findall(gcf,'Type','line');
for i = 1:length(hlines)
    set(hlines(i),'LineWidth',2.0);
    set(hlines(i),'Color','k');
end

if(compartments>1.5)
set(h1,'color',[0 220/255 0] ,'linewidth',4.0);
set(h2,'color',[225/225 102/225 0] ,'linewidth',4.0);
else
    set(h1,'color',[0 0 0] ,'linewidth',4.0);
end
rho_labels = {'0.6' '0.4' '0.2'};
rho_labels2 = {'' '' ''};
ff = findall(f,'type','text');
t=strtrim(get(ff,'String'));
for r=1:length(rho_labels)
    set(ff(strcmp(t,rho_labels{r})),'String',rho_labels2{r})
end

if(any(odour_context>1.5)) %if there's any apple, i.e. if a three-compartment box
    h3 = polar(angles',varargin{3});
    set(h3,'color',[0 0 220/255], 'linewidth', 4.0);
end

%% remove text objects except the main line
delete(findall(ancestor(f,'figure'),'HandleVisibility','off','type','text'));
set(gca,'FontSize',28);
%% save figure
if (strcmpi(save_mode, 'save'))
    set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    if strcmpi(test_type,'s')
        filename = ['_',celltype,'Cell',num2str(cell_number),'Tuning_Test'];
    elseif strcmpi(test_type,'d')
        filename = ['_',celltype,'Cell',num2str(cell_number),'Tuning_Disor'];
    else
        filename = ['_',celltype,'Cell',num2str(cell_number),'Tuning'];
    end
    saveas(gcf,filename, 'fig');
    close(gcf);
end

end

