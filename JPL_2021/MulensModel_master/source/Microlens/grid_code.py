import numpy as np
import scipy.optimize as op
import MulensModel as mm
from MulensModel.utils import Utils
import matplotlib.pyplot as pl
import est_params as dat
import fit_func as fu
class color:
   PURPLE = '\033[1;35;48m'
   CYAN = '\033[1;36;48m'
   BOLD = '\033[1;37;48m'
   BLUE = '\033[1;34;48m'
   GREEN = '\033[1;32;48m'
   YELLOW = '\033[1;33;48m'
   RED = '\033[1;31;48m'
   BLACK = '\033[1;30;48m'
   UNDERLINE = '\033[4;37;48m'
   END = '\033[1;37;0m'



s_step = 70
q_step = 70
# grid_log_s = np.hstack(
#     (np.arange(
#         np.log10(s_minus) - 0.1, np.log10(s_minus) + 0.1, delta_log_s),
#     np.arange(
#         np.log10(s_plus) - 0.1, np.log10(s_plus) + 0.1, delta_log_s)))
# grid_log_s = np.arange(np.log10(.6),np.log10(1.666),delta_log_s)
# grid_log_q = np.arange(-5, -1, delta_log_q)
# print(grid_log_s)
# print(grid_log_q)

grid_log_s = np.linspace(np.log10(.6),np.log10(1.666),s_step)
grid_log_q = np.linspace(-4, -.3,q_step)

grid_save = []
grid = np.empty((5, len(grid_log_s) * len(grid_log_q)))
counter = 0
print('{0:>12} {1:>6} {2:>7} {3:>7} {4:>7}'.format('chi2', 's', 'q', 'alpha', 'rho'))
def grid_fit(tnot,unot,tEnot,anot):
    for log_s in grid_log_s:
        for log_q in grid_log_q:
            # The major and minor images are on opposite sides of the lens:
            if log_s < 0.:
                alpha = anot + 180.0
            else:
                alpha = anot

            # Define the Model and Event
            planet_model = mm.Model({
                't_0': tnot,
                'u_0': unot,
                't_E': tEnot,
                'rho': 10**(-2),
                's': 10.**log_s,
                'q': 10.**log_q,
                'alpha': anot})
            planet_model.set_magnification_methods(dat.magnification_methods)
            planet_event = mm.Event(datasets = [dat.H_data , dat.K_data], model=planet_model)

            # Fit the Event
            result = fu.fit_model(planet_event, parameters_to_fit=['rho', 'alpha'])
            if result.success:
                chi2 = planet_event.get_chi2()
            else:
                chi2 = 40000#np.inf

            # Print and store result of fit
            print('{0:12.2f} {1:6.4f} {2:7.5f} {3:7.2f} {4:7.5f}'.format(
                chi2, 10.**log_s, 10.**log_q,
                planet_event.model.parameters.alpha, planet_event.model.parameters.rho))
            global counter
            grid[0, counter] = log_s
            grid[1, counter] = log_q
            grid[2, counter] = chi2
            grid[3, counter] = planet_event.model.parameters.alpha.value
            grid[4, counter] = planet_event.model.parameters.rho
            counter += 1
            # np.savetxt('test_grid.txt', np.column_stack([grid_save.append(chi2),grid_save.append(log_s),
            # grid_save.append(log_q),grid_save.append(planet_event.model.parameters.alpha.value),
            # grid_save.append(planet_event.model.parameters.rho)]))
    save_chi = np.array(grid[2,:])
    save_s = np.array(grid[0,:])
    save_q = np.array(grid[1,:])
    save_alpha = np.array(grid[3,:])
    save_rho = np.array(grid[4,:])
    index_best = np.argmin(np.array(grid[2,:]))
    index_sorted = np.argsort(np.array(grid[2,:]))
    np.savetxt('grid_70_steps.txt', np.column_stack([save_chi,save_s,save_q,save_alpha,save_rho]),delimiter=' ',fmt='%f')



    n_best = 1
    colors = ['magenta', 'green', 'cyan','yellow','blue', 'red','orange']
    if len(colors) < n_best - 1:
        raise ValueError('colors must have at least n_best -1 entries.')

# Plot the grid

    fig, axes = pl.subplots(nrows=1, ncols=2,figsize=(8,8))
    n_plot = 0
    for i in np.arange(2):
        if i == 0:
            index_logs = np.where(grid_log_s < 0.)[0]
            index_grid = np.where(grid[0, :] < 0.)[0]
        else:
            index_logs = np.where(grid_log_s >= 0.)[0]
            index_grid = np.where(grid[0, :] >= 0.)[0]

        # Plot chi2 map

        chi2 = np.transpose(
                grid[2, index_grid].reshape(len(index_logs), len(grid_log_q)))
        print(chi2)
        im = axes[i].imshow(
            chi2, aspect='auto', origin='lower',
            extent=[
                np.min(grid_log_s[index_logs]),
                np.max(grid_log_s[index_logs]),
                np.min(grid_log_q),
                np.max(grid_log_q)],
            cmap='gray',
            vmin=np.min(grid[2,:]), vmax=np.nanmax(grid[2,np.isfinite(grid[2,:])]))

        # Mark best values: best="X", other good="o"
        if index_best in index_grid:
            axes[i].scatter(grid[0, index_best], grid[1, index_best], marker='x', color='white')
        for j, index in enumerate(index_sorted[1:n_best]):
            if index in index_grid:
                axes[i].scatter(grid[0, index], grid[1, index], marker='o', color=colors[j - 1])

    print(np.min(grid_log_s[index_logs]),np.max(grid_log_s[index_logs]),np.min(grid_log_q),
    np.max(grid_log_q))

    fig.subplots_adjust(right=0.9)
    cbar_ax = fig.add_axes([0.95, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax)

    fig.text(0.5, 0.92, r'$\chi^2$ Map', ha='center')
    fig.text(0.5, 0.04, 'log s', ha='center')
    fig.text(0.04, 0.5, 'log q', va='center', rotation='vertical')
    fig.text(1.1, 0.5, r'$\chi^2$', va='center', rotation='vertical')
    #pl.savefig('grid_1.png', format='png', dpi=500)
    pl.show()

    def make_grid_model(index):
        model = mm.Model({
            't_0': planet_event.model.parameters.t_0,
            'u_0': planet_event.model.parameters.u_0,
            't_E': planet_event.model.parameters.t_E,
            'rho': grid[4, index],
            's': 10.**grid[0, index],
            'q': 10.**grid[1, index],
            'alpha': grid[3, index]})
        model.set_magnification_methods(dat.magnification_methods)
        model_event = mm.Event(datasets = [dat.H_data,dat.K_data], model=model)
        return model

    best_fit_model = make_grid_model(index_best)
    best_fit_event = mm.Event(datasets=[dat.H_data,dat.K_data], model=best_fit_model)
    print(best_fit_event.model)
    print('chi2: {0}'.format(best_fit_event.get_chi2()))
    for j, index in enumerate(index_sorted[1:n_best]):
        model = make_grid_model(index)
        model.plot_lc(t_range= [2457900,2458800],subtract_2450000=True, color=colors[j - 1], lw=2)
        best_fit_event.plot_data(subtract_2450000=True, s=10, zorder=0,marker_list = 'o',markerfacecolor='none')
        pl.xlim(8620,8685)
        print('Other best fit models')
        print(model)

    print(color.RED +'Refining the n_best minima to get the best-fit solution'+ color.END)
    parameters_to_fit = ['t_0', 'u_0', 't_E', 'rho', 'alpha', 's', 'q'] #mcmc

    fits = []
    for index in index_sorted[:n_best]:
        model = make_grid_model(index)
        event = mm.Event(datasets=[dat.H_data,dat.K_data], model=model)
        print(event.model)
        result = fu.fit_model(
        event, parameters_to_fit=parameters_to_fit)
        fits.append([result.fun, result.x])
        print(result)

    chi2 = [x[0] for x in fits]
    fit_parameters = [x[1] for x in fits]
    index_best = np.argmin(chi2)

    # Setup the model and event
    parameters = {}
    for i, parameter in enumerate(parameters_to_fit):
        parameters[parameter] = fit_parameters[index_best][i]
    final_model = mm.Model(parameters)
    final_model.set_magnification_methods(dat.magnification_methods)
    final_event = mm.Event(datasets=[dat.H_data,dat.K_data], model=final_model)
    print(color.RED +'Final best fit model'+ color.END)
    print(final_event.model)
    print('chi2: {0}'.format(final_event.get_chi2()))
    pl.figure(figsize=(10,10))
    pl.axes ([0.09 , 0.08 , 0.9 , 0.9])
    final_event.plot_data (subtract_2450000 = True , label_list =['H_data','K_data'] ,color_list=['red','black'] ,
    zorder_list =[2 , 1] , s =6)
    final_event.plot_model(subtract_2450000 = True , lw=2, color = 'blue',label = 'best fit')
    pl.xlim(8620,8685)
    pl.legend ( loc ='best')
    pl.axes([0.17 , 0.7 , 0.3 , 0.2]) # Figure inset stars here .
    final_model.plot_trajectory (caustics = True)
    pl.xlim ( -0.20 , 0.9)
    pl.ylim ( -0.2 , 0.2)
    pl.show()

    # pl.figure(figsize=(10,10))
    # best_fit_event.plot_model(t_range = [2457900,2458800],subtract_2450000=True, color='black',label='best',zorder=10)
    # best_fit_event.plot_data(subtract_2450000=True, s=5, zorder=10)
    # pl.xlim(8600,8800)
    # #pl.title(filename + ' grid model zoom')
    # #pl.savefig(filename + '_grid_model_zoom', dpi=100)
    # pl.show()
    #
    # pl.figure(figsize=(10,10))
    # best_fit_event.plot_model(t_range = [2457900,2458800],subtract_2450000=True, color='black',label='best',zorder=10)
    # best_fit_event.plot_data(subtract_2450000=True, s=5, zorder=10)
    # pl.xlim(8620,8690)
    #     #pl.title(filename + ' grid model zoom')
    #     #pl.savefig(filename + '_grid_model_zoom', dpi=100)
    # pl.show()
