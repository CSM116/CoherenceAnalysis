function [] = bx_connector(xpos,ypos,color)
    if nargin<3
      color = 'k';
    end
    plot([xpos(1),xpos(1)],[ypos(1)-0.01,ypos(1)],color);
    plot([xpos(1),xpos(2)],[ypos(1),ypos(2)],color);
    plot([xpos(2),xpos(2)],[ypos(2)-0.01,ypos(2)],color);
end