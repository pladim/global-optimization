from ortools.linear_solver import pywraplp

def solve(x21_1, x21_2, x31_1, x31_2, x41_1, x41_2, x32_1, x32_2, x42_1, x42_2, x43_1, x43_2, fa, fv, fu, fb):
    norm21 = x21_1 * x21_1 + x21_2 * x21_2
    norm31 = x31_1 * x31_1 + x31_2 * x31_2
    norm41 = x41_1 * x41_1 + x41_2 * x41_2
    norm32 = x32_1 * x32_1 + x32_2 * x32_2
    norm42 = x42_1 * x42_1 + x42_2 * x42_2
    norm43 = x43_1 * x43_1 + x43_2 * x43_2

    solver = pywraplp.Solver.CreateSolver('GLOP')

    L   = solver.NumVar(0,                  solver.infinity(), 'L')
    ga1 = solver.NumVar(-solver.infinity(), solver.infinity(), 'ga1')
    ga2 = solver.NumVar(-solver.infinity(), solver.infinity(), 'ga2')
    gv1 = solver.NumVar(-solver.infinity(), solver.infinity(), 'gv1')
    gv2 = solver.NumVar(-solver.infinity(), solver.infinity(), 'gv2')
    gu1 = solver.NumVar(-solver.infinity(), solver.infinity(), 'gu1')
    gu2 = solver.NumVar(-solver.infinity(), solver.infinity(), 'gu2')
    gb1 = solver.NumVar(-solver.infinity(), solver.infinity(), 'gb1')
    gb2 = solver.NumVar(-solver.infinity(), solver.infinity(), 'gb2')

    solver.Add(-L * norm21 + (gv1 * x21_1 + gv2 * x21_2) - (ga1 * x21_1 + ga2 * x21_2) <= 0.0)
    solver.Add(-L * norm21 - (gv1 * x21_1 + gv2 * x21_2) + (ga1 * x21_1 + ga2 * x21_2) <= 0.0)
    solver.Add(-0.5 * L * norm21 - (ga1 * x21_1 + ga2 * x21_2) <= fa - fv)
    solver.Add(-0.5 * L * norm21 + (ga1 * x21_1 + ga2 * x21_2) <= fv - fa)
    solver.Add(-0.5 * L * norm21 - (gv1 * x21_1 + gv2 * x21_2) <= fa - fv)
    solver.Add(-0.5 * L * norm21 + (gv1 * x21_1 + gv2 * x21_2) <= fv - fa)

    solver.Add(-L * norm31 + (gu1 * x31_1 + gu2 * x31_2) - (ga1 * x31_1 + ga2 * x31_2) <= 0.0)
    solver.Add(-L * norm31 - (gu1 * x31_1 + gu2 * x31_2) + (ga1 * x31_1 + ga2 * x31_2) <= 0.0)
    solver.Add(-0.5 * L * norm31 - (ga1 * x31_1 + ga2 * x31_2) <= fa - fu)
    solver.Add(-0.5 * L * norm31 + (ga1 * x31_1 + ga2 * x31_2) <= fu - fa)
    solver.Add(-0.5 * L * norm31 - (gu1 * x31_1 + gu2 * x31_2) <= fa - fu)
    solver.Add(-0.5 * L * norm31 + (gu1 * x31_1 + gu2 * x31_2) <= fu - fa)

    solver.Add(-L * norm41 + (gb1 * x41_1 + gb2 * x41_2) - (ga1 * x41_1 + ga2 * x41_2) <= 0.0)
    solver.Add(-L * norm41 - (gb1 * x41_1 + gb2 * x41_2) + (ga1 * x41_1 + ga2 * x41_2) <= 0.0)
    solver.Add(-0.5 * L * norm41 - (ga1 * x41_1 + ga2 * x41_2) <= fa - fb)
    solver.Add(-0.5 * L * norm41 + (ga1 * x41_1 + ga2 * x41_2) <= fb - fa)
    solver.Add(-0.5 * L * norm41 - (gb1 * x41_1 + gb2 * x41_2) <= fa - fb)
    solver.Add(-0.5 * L * norm41 + (gb1 * x41_1 + gb2 * x41_2) <= fb - fa)

    solver.Add(-L * norm32 + (gu1 * x32_1 + gu2 * x32_2) - (gv1 * x32_1 + gv2 * x32_2) <= 0.0)
    solver.Add(-L * norm32 - (gu1 * x32_1 + gu2 * x32_2) + (gv1 * x32_1 + gv2 * x32_2) <= 0.0)
    solver.Add(-0.5 * L * norm32 - (gv1 * x32_1 + gv2 * x32_2) <= fv - fu)
    solver.Add(-0.5 * L * norm32 + (gv1 * x32_1 + gv2 * x32_2) <= fu - fv)
    solver.Add(-0.5 * L * norm32 - (gu1 * x32_1 + gu2 * x32_2) <= fv - fu)
    solver.Add(-0.5 * L * norm32 + (gu1 * x32_1 + gu2 * x32_2) <= fu - fv)

    solver.Add(-L * norm42 + (gb1 * x42_1 + gb2 * x42_2) - (gv1 * x42_1 + gv2 * x42_2) <= 0.0)
    solver.Add(-L * norm42 - (gb1 * x42_1 + gb2 * x42_2) + (gv1 * x42_1 + gv2 * x42_2) <= 0.0)
    solver.Add(-0.5 * L * norm42 - (gv1 * x42_1 + gv2 * x42_2) <= fv - fb)
    solver.Add(-0.5 * L * norm42 + (gv1 * x42_1 + gv2 * x42_2) <= fb - fv)
    solver.Add(-0.5 * L * norm42 - (gb1 * x42_1 + gb2 * x42_2) <= fv - fb)
    solver.Add(-0.5 * L * norm42 + (gb1 * x42_1 + gb2 * x42_2) <= fb - fv)

    solver.Add(-L * norm43 + (gb1 * x43_1 + gb2 * x43_2) - (gu1 * x43_1 + gu2 * x43_2) <= 0.0)
    solver.Add(-L * norm43 - (gb1 * x43_1 + gb2 * x43_2) + (gu1 * x43_1 + gu2 * x43_2) <= 0.0)
    solver.Add(-0.5 * L * norm43 - (gu1 * x43_1 + gu2 * x43_2) <= fu - fb)
    solver.Add(-0.5 * L * norm43 + (gu1 * x43_1 + gu2 * x43_2) <= fb - fu)
    solver.Add(-0.5 * L * norm43 - (gb1 * x43_1 + gb2 * x43_2) <= fu - fb)
    solver.Add(-0.5 * L * norm43 + (gb1 * x43_1 + gb2 * x43_2) <= fb - fu)

    solver.Minimize(L)
    status = solver.Solve()
    if status == pywraplp.Solver.OPTIMAL:
        print('OPTIMAL')
    else:
        print('NOT OPTIMAL')

    return solver.Objective().Value()
