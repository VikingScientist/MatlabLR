function lrgui
% LRGUI Select a data set from the pop-up menu, then
% click one of the plot-type push buttons. Clicking the button
% plots the selected data in the axes.



  %  Create and then hide the GUI as it is being constructed.
  f = figure('Visible','off','Position',[360,500,450,285]);

  %  Construct the components.
  hnew = uicontrol('Style','pushbutton','String','New',...
         'Position',[315,240,70,25],...
         'Callback',@newbutton_Callback);
  hrefine = uicontrol('Style','pushbutton',...
         'String','Refine',...
         'Position',[315,200,70,25],...
         'Callback',@refinebutton_Callback); 
  hload = uicontrol('Style','pushbutton','String','Load',...
         'Position',[315,160,70,25],...
         'Callback',@loadbutton_Callback);
  hsave = uicontrol('Style','pushbutton','String','Save',...
         'Position',[315,120,70,25],...
         'Callback',@savebutton_Callback);
  htext = uicontrol('Style','text','String','Refinement strategy',...
         'Position',[300,60,80,35]);
  hpopup = uicontrol('Style','popupmenu',...
         'String',{'Elements','Functions'},...
         'Position',[300,30,100,25],...
         'Callback',@popup_menu_Callback);
  ha = axes('Units','Pixels','Position',[50,60,200,185]); 
  align([hnew,hrefine,hload,hsave,htext,hpopup],'Center','None');
  
  % Create the data to plot.
  current_mesh = LRSplineSurface([3,3], [10,10]);
  
  % Initialize the GUI.
  % Change units to normalized so components resize 
  % automatically.
  f.Units = 'normalized';
  ha.Units = 'normalized';
  hnew.Units = 'normalized';
  hload.Units = 'normalized';
  hsave.Units = 'normalized';
  hrefine.Units = 'normalized';
  htext.Units = 'normalized';
  hpopup.Units = 'normalized';

  % add key-listener for keyboard hotkeys
  set(gcf, 'KeyPressFcn', @keypress);
  
  %Create a mesh in the axes.
  H = current_mesh.plot('parametric');
  ref_elements = true;
  refine_list = [];

  redraw_mesh();

  % assign the gui a name to appear in the window title.
  f.Name = 'lr mesh generation';
  % move the gui to the center of the screen.
  movegui(f,'center')
  % make the gui visible.
  f.Visible = 'on';

  %  Callbacks for lrgui. These callbacks automatically
  %  have access to component handles and initialized data 
  %  because they are nested at a lower level.

  %  Pop-up menu callback. Read the pop-up menu Value property
  %  to determine which item is currently displayed and make it
  %  the current data.
    function popup_menu_Callback(source,eventdata)
      % Determine the selected data set.
      str = source.String;
      val = source.Value;
      % Set current data to the selected data set.
      switch str{val};
      case 'Elements' % User selects Elements.
        ref_elements = true;
      case 'Functions' % User selects Functions.
        ref_elements = false;
      end
      redraw_mesh();
    end
 
  % Push button callbacks. 

  function newbutton_Callback(source,eventdata) 
  % Display surf plot of the currently selected data.
    % surf(current_data);
    x = inputdlg({'p1', 'p2', 'n1', 'n2'}, 'New LR spline', [1 8], {'3' '3' '10' '10'});
    if numel(x) > 0
      p = [str2num(x{1}), str2num(x{2})];
      n = [str2num(x{3}), str2num(x{4})];
      current_mesh = LRSplineSurface(p,n);
      redraw_mesh();
    end
  end

  function refinebutton_Callback(source,eventdata) 
  % Display contour plot of the currently selected data.
    % contour(current_data);
    if ref_elements
      current_mesh.refine(refine_list, 'elements');
    else
      current_mesh.refine(refine_list, 'basis');
    end
    redraw_mesh();
  end 

  function loadbutton_Callback(source,eventdata) 
  % Display mesh plot of the currently selected data.
    [file path] = uigetfile({'*.g2;*.lr', 'Geometry files (*.g2,*.lr)'});
    if file ~= 0
      current_mesh.load([path,file]);
      redraw_mesh();
    end
  end

  function savebutton_Callback(source,eventdata) 
  % Display mesh plot of the currently selected data.
    [file path] = uiputfile('out.lr');
    if file ~= 0
      current_mesh.save([path,file]);
    end
  end

  % Keyboard callbacks (hotkeys)

  function keypress(src, e)
    switch e.Key
      case 'r'
        refinebutton_Callback(src, e);
      case 'e'
        hpopup.Value = 1;
        popup_menu_Callback(hpopup,e);
      case 'f'
        hpopup.Value = 2;
        popup_menu_Callback(hpopup,e);
      case 'l'
        loadbutton_Callback(src, e);
      case 's'
        savebutton_Callback(src, e);
      case 'n'
        newbutton_Callback(src, e);
    end
  end

  % Mouse callbacks

  function mouseclick(src,e)
    x = e.IntersectionPoint;
    u = x(1);
    v = x(2);
    if ref_elements
      i  = current_mesh.getElementContaining(u,v);
      u0 = current_mesh.elements(i,1);
      v0 = current_mesh.elements(i,2);
      u1 = current_mesh.elements(i,3);
      v1 = current_mesh.elements(i,4);
      hold on;
      patch([u0,u1,u1,u0], [v0,v0,v1,v1], [.8, .8, .8]);
      refine_list = [refine_list, i];
    else
      i  = pick_basisfunction(u,v);
      n = size(current_mesh.knots,1);
      dots = get(gca,'children');
      dots(i).MarkerFaceColor = [1,0,0];
      refine_list = [refine_list, n+1-i];
    end
  end

  function i = pick_basisfunction(u,v)
    mindist = inf;
    j = 1;
    for l=get(gca, 'children')'
      if numel(l.MarkerFaceColor) == 3 % this is (rgb)-values or 'none'-string
        dist = norm([u,v] - [l.XData, l.YData]);
        if dist < mindist
          i = j;
          mindist = dist;
        end
      end
      j = j+1;
    end 
  end

  function redraw_mesh()
    refine_list = [];
    cla
    hold on;
    if ref_elements
      current_mesh.plot('parametric');
    else
      current_mesh.plot('parametric', 'basis');
    end
    for line = get(gca, 'children')'
      set(line, 'HitTest', 'off');
    end
    set(gca, 'ButtonDownFcn', @mouseclick);
  end

end 
