function [] = bx_connector(xpos,ypos,color)
    if nargin<3
      color = 'w';
    end
%     style = ':';
%     linew = 1.175;
    style = '-';
    linew = 0.6;
    plot([xpos(1),xpos(1)],[ypos(1)-0.005,ypos(1)],strcat(style,color),'LineWidth',linew);
    plot([xpos(1),xpos(2)],[ypos(1),ypos(2)],strcat(style,color),'LineWidth',linew);
    plot([xpos(2),xpos(2)],[ypos(2)-0.005,ypos(2)],strcat(style,color),'LineWidth',linew);
end