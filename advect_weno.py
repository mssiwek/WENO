import numpy as np
import matplotlib.pyplot as plt
from weno import *
from test_functions import *

wavespeed        = 1.0
cfl_number       = 0.5
func_id          = (input('Which IC function do you want to advect? Examples include: gaussian, simple_step, square_well, quadratic_well, trig_disc. See test_functions.py for more info. \n'))
func             = eval(func_id)
num_cells        = eval(input("At which resolution? \n"))
plot_L2          = input('Do you want to plot the L2 error? This may take some time. [y/n]\n')
if plot_L2 == 'y':
    all_res = eval(input('For which resolutions do you want to plot L2? (give in format: [res0,res1,res3,...]) \n'))

t_final         = eval(input('End time:\n'))
first_order      = input("Do you also want to plot the first order solution for comparison? [y/n]\n")

def solution_at_time(t, x, func, a = wavespeed):
    fmap             = lambda f: lambda x: [f(xi) for xi in x]
    wrap             = lambda x: wrap(x + 2) if x < -1 else (wrap(x - 2) if x > 1 else x)
    compose          = lambda f, g: lambda x: f(g(x))
    solution_at_time = lambda t: lambda x: compose(func, wrap)(x - a*t)
    return(np.array(fmap(solution_at_time(t))(x)))
    #return 1.5 + np.sin(2 * np.pi * (xc - wavespeed * t))

def compute_du_upwind(xv, u, dt):
    dx = np.diff(xv)
    f0 = wavespeed * u
    fhat = f0[:-1] if wavespeed > 0.0 else f0[1:]
    return -np.diff(fhat) * dt / dx


def compute_du_weno(xv, u, dt):
    dx = np.diff(xv)
    f0 = wavespeed * u

    def fhat_at_iph(i):
        if wavespeed > 0.0:
            fi = f0[i-3:i+2]
            return nonlinear_weighted(*fi)
        else:
            assert(false)

    fhat = np.array([fhat_at_iph(i) for i in range(3, len(u) - 2)])
    return -(fhat[1:] - fhat[:-1]) * dt / dx



def extend(u, num_guard_zones):
    return np.concatenate([u[-num_guard_zones:], u, u[0:num_guard_zones]])



def update(t, u, xv, compute_du = compute_du_weno):
    if compute_du == compute_du_weno:
        ng = 3
    if compute_du == compute_du_upwind:
        ng = 1

    dx = np.diff(xv)
    dt = abs(cfl_number * dx[0] / wavespeed)

    k1 = compute_du(xv, extend(u +      0.0, ng), dt)
    k2 = compute_du(xv, extend(u + k1 * 0.5, ng), dt)
    k3 = compute_du(xv, extend(u + k2 * 0.5, ng), dt)
    k4 = compute_du(xv, extend(u + k3 * 1.0, ng), dt)
    t += dt
    u += (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return t, u



def evolve_with_num_cells(num_cells, func, compute_du = compute_du_weno):
    print(num_cells)
    xv = np.linspace(-1, 1, num_cells + 1)
    dx = np.diff(xv)
    xc = 0.5 * (xv[1:] + xv[:-1])
    u = solution_at_time(0.0, xc, func)
    #print("u:", u)
    t = 0.0

    us = u
    ts = [t]

    #define tnext here
    #tnext = 0
    while t < t_final:
        t, u = update(t, u, xv, compute_du = compute_du)
        if (t-ts[-1]) >= 0.01:
            #print("t:", t)
            us = np.vstack((us, u))
            ts.append(t)
            #tnext += 0.001

    L2 = np.sum(((u - solution_at_time(t, xc, func))**2)*dx)**0.5
    print("At resolution = {}, L2 error = {:04e}".format(num_cells, L2))
    return ts, xc, us, L2



def plot_solution(ts, xs, us, ts_1 = None, xs_1 = None, us_1 = None):
    n = 1
    for j in range(1,len(ts)):
        # Plotting settings
        import matplotlib as mpl
        mpl.rc('font', **{'family': 'serif', 'sans-serif': ['Times']})
        mpl.rc('lines', solid_capstyle='round')
        mpl.rc('mathtext', fontset='cm')
        plt.rcParams.update({'grid.alpha': 0.5})
        mpl.rcParams['font.size'] = 50
        if first_order == 'n':
            fig = plt.figure(figsize=(30,15))
            ax1 = fig.add_subplot(111)
            ax1.plot(xs, solution_at_time(ts[j], xs, func), linewidth = 7, label = "Analytic, t=%.2f" %ts[j], zorder=1)
            ax1.plot(xs, us[j], '.', markersize='50', label = "Numeric, t=%.2f" %ts[j], zorder = 0)
            ax1.set_title(r'$\rm N_{cells} = %04d, \ with \ WENO$' %len(xs), y=1.05)
            ax1.legend(loc = 'lower right')
            ax1.set_xlabel(r'$\rm x$', fontsize = 60)
            ax1.set_ylabel(r'$\rm \rho$', labelpad = 35, fontsize = 60)
            plt.tight_layout()
            plt.savefig('plots/' + func_id+'_%03d' %j)
            plt.close()
        if first_order == 'y':
            fig = plt.figure(figsize=(30,30))
            ax = fig.add_subplot(111)
            ax.spines['top'].set_color('none')
            ax.spines['bottom'].set_color('none')
            ax.spines['left'].set_color('none')
            ax.spines['right'].set_color('none')
            ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

            ax1 = fig.add_subplot(211)
            ax2 = fig.add_subplot(212)
            ax2.plot(xs_1, solution_at_time(ts_1[j], xs_1, func), linewidth = 7, label = "Analytic, t=%.2f" %ts_1[j], zorder=1)
            ax2.plot(xs_1, us_1[j], '.', markersize='50', label = "Numeric, t=%.2f" %ts[j], zorder = 0)
            ax2.set_title(r'$\rm N_{cells} = %04d, \ without \ WENO$' %len(xs), y=1.05)
            ax2.legend(loc = 'lower right')
            #ax2.set_xlabel(r'$\rm x$', fontsize = 60)
            #ax2.set_ylabel(r'$\rm \rho$', labelpad = 35, fontsize = 60)

            ax1.plot(xs, solution_at_time(ts[j], xs, func), linewidth = 7, label = "Analytic, t=%.2f" %ts[j], zorder=1)
            ax1.plot(xs, us[j], '.', markersize='50', label = "Numeric, t=%.2f" %ts[j], zorder = 0)
            ax1.set_title(r'$\rm N_{cells} = %04d, \ with \ WENO$' %len(xs), y=1.05)
            ax1.legend(loc = 'lower right')
            #ax1.set_xlabel(r'$\rm x$', fontsize = 60)
            #ax1.set_ylabel(r'$\rm \rho$', labelpad = 35, fontsize = 60)

            ax.set_xlabel(r'$\rm x$', fontsize = 60)
            ax.set_ylabel(r'$\rm \rho$', labelpad = 35, fontsize = 60)
            plt.tight_layout()
            plt.savefig('plots/' + func_id+'_%03d' %j)
            plt.close()

    import os
    cmd = 'ffmpeg -r 10 -f image2 -s 1920x1080 -i plots/%s_' %func_id + r'%03d.png' + ' -vcodec libx264 -crf 25  -pix_fmt yuv420p videos/%s.mp4' %func_id
    os.system(cmd)



def L2_plot(all_n_cells, L2s):
    import matplotlib as mpl
    mpl.rc('font', **{'family': 'serif', 'sans-serif': ['Times']})
    mpl.rc('lines', solid_capstyle='round')
    mpl.rc('mathtext', fontset='cm')
    plt.rcParams.update({'grid.alpha': 0.5})
    mpl.rcParams['font.size'] = 50

    fig = plt.figure(figsize=(30,20))
    ax1 = fig.add_subplot(111)
    lw = 7
    ax1.loglog(all_n_cells, L2s, linewidth = lw, label = r'$\rm %s \ convergence$' %func_id)
    ax1.loglog(all_n_cells, (np.array(all_n_cells))**(-5.)*(L2s[1]*all_n_cells[1]**5), linewidth = lw, label = r'$\rm L_2 \propto \delta x ^5$')
    ax1.set_xlabel(r'$\rm N_{cells}$')
    ax1.set_ylabel(r'$\rm L_{2}$')
    ax1.legend()
    plt.savefig('plots/' + func_id+ 'L2.png')
    plt.close()


if __name__ == "__main__":

    if plot_L2 == 'y':
        L2s = []
        for n_cells in all_res:
            print('n_cells', n_cells)
            ts, xs, us, L2 = evolve_with_num_cells(n_cells, func, compute_du = compute_du_weno)
            L2s.append(L2)
        np.savetxt('L2.txt', np.array(L2s))
        L2_plot(all_res,L2s)

    if first_order == 'n':
        ts, xs, us, L2 = evolve_with_num_cells(num_cells, func, compute_du = compute_du_weno)
        plot_solution(ts, xs, us)
    if first_order == 'y':
        ts, xs, us, L2 = evolve_with_num_cells(num_cells, func, compute_du = compute_du_weno)
        ts_1, xs_1, us_1, L2_1 = evolve_with_num_cells(num_cells, func, compute_du = compute_du_upwind)
        plot_solution(ts, xs, us, ts_1 = ts_1, xs_1 = xs_1, us_1 = us_1)
