function u = sapprox(u_00, A_mas, G, dt, t_f, tol, max_iter)
    u = u_00;
    prev_u = u_00;

    while t <= t_f
        error_val = 10;
        u_a = u;
        p = 0;

        while error_val < tol && p < max_iter
            A = A_mas(u);
            v = G(u);
            u = A \ (u_k + v);
            error_val = norm(u - u_a, inf) / norm(u, inf);
            u_a = u;
            p += 1;
        endwhile

        u_k = u;
        t += dt;
    endwhile
endfunction