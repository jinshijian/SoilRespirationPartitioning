function plot_global_map(lats, lons, rc_total)
% sey projection and grid
axis equal
m_proj('robinson','lat',[-62 90],'lon',[-180 180], 'linewidth', 1.5); % robinson Mollweide
hold on
m_grid('tickdir','in','linestyle','none','backcolor',[.9 .99 1], 'xticklabels',[], ...
    'fontsize',8);

% plot image
im = m_image(lons,lats,rc_total);

% add border
M=m_shaperead('landareas'); 
for k=1:length(M.ncst)    
     m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'color','k'); 
end 
shading flat;
view(0,90);
hold off