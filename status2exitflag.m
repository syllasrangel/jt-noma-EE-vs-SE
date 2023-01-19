function EXIT_FLAG = status2exitflag(status)
    switch(status)
        case 'Solved'
            EXIT_FLAG = 1;
        case 'Suboptimal'
            EXIT_FLAG = 2;
        case 'Inaccurate/Solved'
            EXIT_FLAG = 3;
        case 'Infeasible'
            EXIT_FLAG = -1;
        case 'Failed'
            EXIT_FLAG = -2;
        case 'Unbounded'
            EXIT_FLAG = -3;
        case 'Overdetermined'
            EXIT_FLAG = -4;
        otherwise
            EXIT_FLAG = 0;
    end
end