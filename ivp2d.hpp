#ifndef IVP2D_H_

#define IVP2D_H_

#include <functional>
#include "octave/oct.h"
#include "octave/parse.h"

#define GMRES(x, A, b, k, tol, max_iter) \
    octave_value_list GMRESinputList; \
    GMRESinputList(0) = A; \
    GMRESinputList(1) = b; \
    GMRESinputList(2) = k; \
    GMRESinputList(3) = tol; \
    GMRESinputList(4) = max_iter; \
    x = (octave::feval("gmres", GMRESinputList))

#define NORM(n, v) \
    octave_value_list normList; \
    normList(0) = v;
    normList(1) = "inf";
    n = (octave::feval("norm", normList))

ColumnVector explicit_ivp(ColumnVector u_0, std::function G, float dt, float tf)
{
    ColumnVector u(u_0);
    float t = dt;

    while (t <= tf)
    {
        u = u + G(u);
        t += dt;
    }

    return u;
}

ColumnVector sapprox(ColumnVector u_0, std::function A_mas, std::function G, float dt, float tf, float tol, int max_iter)
{
    ColumnVector u(u_0);
    ColumnVector prev_u(u_0);
    float t = dt;

    while (t <= tf)
    {
        float error = 10.0f;
        ColumnVector u_a(u);
        
        for (int i = 0; i < max_iter && error < tol; ++i)
        {
            Matrix A = A_mas(u);
            ColumnVector v = G(u);
            GMRES(u, A, prev_u + v, 200, tol, max_iter);
            float a = 0.0f, b = 0.0f;
            NORM(a, u - u_a);
            NORM(b, u);
            error = a / b;
            u_a = u;
        }

        prev_u = u;
        t += dt;
    }

    return u;
}

ColumnVector newton(ColumnVector u_0, std::function A, std::function F, float dt, float tf, float tol, int max_iter)
{
    ColumnVector u(u_0);
    float t = dt;

    while (t <= tf)
    {
        ColumnVector f = F(u);
        float r;
        NORM(r, f);
        float r_0 = r;
        float delta = tol * r_0;

        for (int i = 0; i < max_iter && r > delta; ++i)
        {
            Matrix a = A(u);
            f = F(u);
            ColumnVector s;
            GMRES(s, a, -f, 200, tol, max_iter);
            u = u + s;
            NORM(r, f);
        }

        t += dt;
    }

    return u;
}

#endif