function u = explicit(u_0, G, dt, t_f)
    u = u_0;
    t = dt;
    while t < t_f
        u += G(u);
        t += dt;
    endwhile
endfunction