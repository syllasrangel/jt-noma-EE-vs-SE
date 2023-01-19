function [BS1_position, BS2_position, inner_users_pos_BS1, inner_users_pos_BS2, JT_users] = generate_scenario_1(N_inner_users, N_JT_users, cell_radius, inner_circle_radius, inner_users_dist)
% Generates a scenario with N_inner_users inner users and N_JT_users edge
% users, where the edge users are randomly distributed in the intersection
% area between the two BSs and the inner users are randomly distributed if
% inner_users_dist is zeros(N_inner_users,1). Otherwise, the inner users
% positions are acordding to the distances from its respective serving BS
% defined in inner_users_dist.
%
%generate_scenario_1(N_inner_users, N_JT_users, cell_radius, inner_circle_radius, inner_users_dist)

%------- CoMP-NOMA Deployment Scenario–1 -------
    BS1_position = 0 + 1i*0;
    BS2_position = (cell_radius + inner_circle_radius) + 1i*0;

%     % inner users with fixed distance
%     %d_u2 = 0.3;
%     d_u2 = 0.125;
%     a = rand*2*pi;
%     b = rand*2*pi;
%     inner_users_pos_BS1 = [d_u1(d_u1_idx)*(cos(a)+1i*sin(a)) + BS1_position; d_u2*(cos(b)+1i*sin(b)) + BS1_position];
%     a = rand*2*pi;
%     b = rand*2*pi;
%     inner_users_pos_BS2 = [d_u1(d_u1_idx)*(cos(a)+1i*sin(a)) + BS2_position; d_u2*(cos(b)+1i*sin(b)) + BS2_position];
    
    
    % Generate two cells with the specified radius and randomly distribued
    % inner and edge users.
    [inner_users_pos_BS1, edge_users_pos_BS1] = generate_cell(BS1_position, N_inner_users, N_JT_users*5, cell_radius, inner_circle_radius);
    [inner_users_pos_BS2, edge_users_pos_BS2] = generate_cell(BS2_position, N_inner_users, N_JT_users*5, cell_radius, inner_circle_radius);
    
    if(any(inner_users_dist) && length(inner_users_dist) == N_inner_users)
        inner_users_pos_BS1 = zeros(N_inner_users,1);
        inner_users_pos_BS2 = zeros(N_inner_users,1);
        for d_idx = 1:length(inner_users_dist)
            a = rand*2*pi;
            b = rand*2*pi;
            inner_users_pos_BS1(d_idx) = inner_users_dist(d_idx)*(cos(a)+1i*sin(a)) + BS1_position;
            inner_users_pos_BS2(d_idx) = inner_users_dist(d_idx)*(cos(b)+1i*sin(b)) + BS2_position;
        end
    elseif (any(inner_users_dist))
        error("Invalid parameter 'inner_users_dist'.");
    end

    % Finds the edge users that can be served by both BSs
    JT_users = zeros(N_JT_users,1);
    while(~all(JT_users))
        users_cell_1 = abs(edge_users_pos_BS1 - BS2_position) <= cell_radius;
        users_cell_2 = abs(edge_users_pos_BS2 - BS1_position) <= cell_radius;
        if(any(users_cell_1))
            users_idx = find(users_cell_1);
            for ii = 1:length(users_idx)
                zero_idx = find(~abs(JT_users),1);
                if(zero_idx)
                    JT_users(zero_idx) = edge_users_pos_BS1(users_idx(ii));
                end
            end
        end
        if(any(users_cell_2))
            users_idx = find(users_cell_2);
            for ii = 1:length(users_idx)
                zero_idx = find(~abs(JT_users),1);
                if(zero_idx)
                    JT_users(zero_idx) = edge_users_pos_BS2(users_idx(ii));
                end
            end
        end
        
        if(~all(JT_users))
            [~, edge_users_pos_BS1] = generate_cell(BS1_position, N_inner_users, N_JT_users*5, cell_radius, inner_circle_radius);
            [~, edge_users_pos_BS2] = generate_cell(BS2_position, N_inner_users, N_JT_users*5, cell_radius, inner_circle_radius);
        end
    end

%     if(sum(users_cell_1)>0)
%         users_idx = find(users_cell_1);
%         JT_CoMP_user = edge_users_pos_BS1(users_idx(1));
%     else
%         users_idx = find(users_cell_2);
%         JT_CoMP_user = edge_users_pos_BS2(users_idx(1));
%     end
end