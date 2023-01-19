function idx = two_dim_2_one_dim(j, bs, N_users, isJT)
    if(j > N_users || (~isJT && bs > 1 && j > N_users-1))
        error("User index out of range")
    end
    if (~isJT)
        idx = j + (N_users)*(bs-1);
        if(bs>2)
            idx = idx - (bs-2);
        end
    else
        idx = j + (N_users)*(bs-1);
    end
end