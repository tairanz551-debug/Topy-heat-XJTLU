function fast_convergence_demo()
    % Create figure window
    fig = figure('Position',[200 200 800 600]);
    
    % Create axes
    ax = axes('Parent',fig,'Position',[0.1 0.3 0.8 0.6]);
    hold(ax,'on'); grid(ax,'on');
    xlabel(ax,'x'); ylabel(ax,'f(x)');
    title(ax,'Fast Convergence Function Comparison');
    
    % Generate x values
    x = linspace(0,1,1000);
    
    % Define function set
    funcs = {
        {@(x,n) x.^n, 'Power function x^n', [0.1 10], 5};
        {@(x,k) (exp(k*x)-1)/(exp(k)-1), 'Exponential function (e^{kx}-1)/(e^k-1)', [1 20], 10};
        {@(x,k) 1-exp(-k*x), 'Modified exponential 1-e^{-kx}', [1 20], 10};
        {@(x,k) exp(-k*(1-x).^2), 'Composite function e^{-k(1-x)^2}', [1 20], 5};
    };
    
    % Color set
    colors = lines(length(funcs));
    
    % Plot initial curves
    plots = gobjects(length(funcs),1);
    for i = 1:length(funcs)
        y = funcs{i}{1}(x,funcs{i}{4});
        plots(i) = plot(ax,x,y,'Color',colors(i,:),'LineWidth',2,...
                       'DisplayName',funcs{i}{2});
    end
    legend(ax,'show','Location','northwest');
    
    % Add control panel
    uicontrol('Style','text','Position',[50 80 120 20],...
              'String','Select function type:');
    popup = uicontrol('Style','popupmenu','Position',[50 60 150 20],...
                     'String',funcs(:,2),'Value',1,...
                     'Callback',@update_controls);
    
    uicontrol('Style','text','Position',[50 30 120 20],...
              'String','Parameter value:');
    slider = uicontrol('Style','slider','Position',[50 10 150 20],...
                      'Min',1,'Max',20,'Value',5,...
                      'Callback',@update_plot);
    txt = uicontrol('Style','text','Position',[210 10 50 20],...
                   'String','5.0');
    
    % Update functions
    function update_controls(~,~)
        idx = popup.Value;
        slider.Enable = 'on';
        range = funcs{idx}{3};
        slider.Min = range(1);
        slider.Max = range(2);
        slider.Value = funcs{idx}{4};
        txt.String = num2str(slider.Value);
        update_plot();
    end
    
    function update_plot(~,~)
        idx = popup.Value;
        param = slider.Value;
        txt.String = num2str(param);
        
        for i = 1:length(funcs)
            if i == idx
                y = funcs{i}{1}(x,param);
                set(plots(i),'YData',y);
                plots(i).Visible = 'on';
            else
                plots(i).Visible = 'off';
            end
        end
    end
end