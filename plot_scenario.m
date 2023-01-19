function plot_scenario(BS_position_1, inner_users_pos_1, edge_users_pos_1, BS_position_2, inner_users_pos_2, edge_users_pos_2, cell_radius, inner_circle_radius)

figure,
plot(inner_users_pos_1,'.b','MarkerSize',20),
hold on,
plot(edge_users_pos_1,'.g','MarkerSize',20),
plot(real(BS_position_1),imag(BS_position_1),'*r'),
text(real(BS_position_1),imag(BS_position_1),'\leftarrow BS 1'),
plot_circle(real(BS_position_1),imag(BS_position_1),cell_radius),
plot_circle(real(BS_position_1),imag(BS_position_1),inner_circle_radius),

plot(inner_users_pos_2,'.b','MarkerSize',20),
plot(edge_users_pos_2,'.y','MarkerSize',20),
plot(real(BS_position_2),imag(BS_position_2),'*r'),
text(real(BS_position_2),imag(BS_position_2),'\leftarrow BS 2'),
plot_circle(real(BS_position_2),imag(BS_position_2),cell_radius),
plot_circle(real(BS_position_2),imag(BS_position_2),inner_circle_radius),
%axis square,
axis equal,

% Finds the users that can be served by both BSs
users_cell_1 = abs(edge_users_pos_1 - BS_position_2) <= cell_radius;
users_cell_2 = abs(edge_users_pos_2 - BS_position_1) <= cell_radius;

if(sum(users_cell_1)>0)
    plot(edge_users_pos_1(users_cell_1),'oc','MarkerSize',10),
end
if(sum(users_cell_2)>0)
    plot(edge_users_pos_2(users_cell_2),'oc','MarkerSize',10),
end
end