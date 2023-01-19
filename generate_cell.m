function [inner_users_pos, edge_users_pos] = generate_cell(cell_center, N_inner_users, N_edge_users, cell_radius, inner_circle_radius)
%Input: Number of inner users, Number of edge users, Cell radius and inner
%circle radius
%
%Output: inner user positions and edge user positions

inner_users_pos = zeros(N_inner_users,1);
edge_users_pos = zeros(N_edge_users,1);

for j=1:N_inner_users
    a = rand*2*pi;
    r = cell_radius*sqrt(rand);
    while(r > inner_circle_radius)
        r = cell_radius*sqrt(rand);
    end
    x = r*cos(a);
    y = r*sin(a);
    inner_users_pos(j) = x+1i*y + cell_center;
end

for j=1:N_edge_users
    a = rand*2*pi;
    r = cell_radius*sqrt(rand);
    while(r < inner_circle_radius)
        r = cell_radius*sqrt(rand);
    end
    x = r*cos(a);
    y = r*sin(a);
    edge_users_pos(j) = x+1i*y + cell_center;
end

end