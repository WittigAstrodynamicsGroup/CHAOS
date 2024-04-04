#python code to define plotting functions to be used in statistical analysis.

'''
This code will automate the creation of plots, using  both matplotlib and plotly

'''
import matplotlib.pyplot as plt
import numpy as np
import scipy
# import plotly.graph_objects as go
from scipy.optimize import curve_fit



'''
Code to define function to create the histograms for the statistical plots of M-C simulations.
The data output can be used with either plotly or matplotlib.pyplot as the user decides. 
'''

def exp(x, *p):
    #exponential function for curve fitting -- A*X**n
    A, n = p
    return A * x**n

def gauss(x, *p):
    '''
    Gaussian distribution function
    '''
    A, mu, sigma = p
    nominator = (x - mu)**2
    denominator = 2 * (sigma)**2
    return A*np.exp(-(nominator / denominator))

def chi(x, k):
    ((x**(k-1)) *  np.exp(-(x**2)/2)) / (2**((k/2) - 1 ) * scipy.special.gamma(k/2))

def hist_ply(data, mean, std_dev, nbins):
    hist, bin_edges = np.histogram(data,  bins = nbins, density=False)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    start = bin_edges[0]
    end = bin_edges[-1]
    # Define model function to be used to fit to the data above:
    
    return(bin_centres, hist, start, end)

def hist_ply_chi(data, mean, std_dev, nbins):
    hist, bin_edges = np.histogram(data,  bins = nbins, density=False)
    bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
    start = bin_edges[0]
    end = bin_edges[-1]
    # Define model function to be used to fit to the data above:
    '''
    p0 is the initial guess for the fitting coefficients (A, mu and sigma above).
    we use the std dev and mean as calculated above.
    '''
    p0 = [3]

    coeff, var_matrix = curve_fit(chi, bin_centres, hist, p0 = p0)

    # Get the fitted curve
    hist_fit = gauss(bin_centres, *coeff)
    return(bin_centres, hist_fit, start, end)

def individual_timetracks(x_values, *data, save_path, save=1, show=0):
    time_x = x_values
    path_save_stat = save_path
    nerrorbars = int(len(time_x) / 10)
    print('total length:', len(time_x))
    print('error every:', nerrorbars)
    omega_x_track, omega_x_err, omega_y_track, omega_y_err, omega_z_track, omega_z_err, omega_norm_track, omega_norm_err = data
    plt.errorbar(time_x, omega_x_track, yerr=omega_x_err, errorevery=nerrorbars, label='mean x-value of angular velocity', capsize=5, ecolor='#D4A368')
    
    plt.legend()
    plt.title('Time profile of the x-component angular velocity with standard deviation')
    plt.xlabel('time [days]')
    plt.ylabel('angular velocity [rad/s]')
    if save == 1:
        plt.savefig(path_save_stat + 'time_track_x.pdf')
    elif show == 1:
        plt.show()
    plt.clf()

    
    plt.errorbar(time_x, omega_y_track, yerr=omega_y_err, errorevery=200, label ='mean y-value of angular velocity', capsize=5, ecolor='#D4A368')
    plt.legend()
    plt.title('Time profile of the y-component angular velocity with standard deviation')
    plt.xlabel('time [days]')
    plt.ylabel('angular velocity [rad/s]')
    if save == 1:
        plt.savefig(path_save_stat + 'time_track_y.pdf')
    elif show == 1:
        plt.show()
    plt.clf()

    
    plt.errorbar(time_x, omega_z_track, yerr=omega_z_err, errorevery=200, label ='mean z-value of angular velocity', capsize=5, ecolor='#D4A368')
    plt.legend()
    plt.title('Time profile of the z-component angular velocity with standard deviation')
    plt.xlabel('time [days]')
    plt.ylabel('angular velocity [rad/s]')
    if save == 1:
        plt.savefig(path_save_stat + 'time_track_z.pdf')
    elif show == 1:
        plt.show()
    plt.clf()


    plt.errorbar(time_x, omega_norm_track, yerr=omega_norm_err, errorevery=200, label='mean angular velocity norm', capsize=5, ecolor='#D4A368')
    plt.legend()
    plt.title('Time profile of the angular velocity norm with standard deviation')
    plt.xlabel('time [days]')
    plt.ylabel('angular velocity [rad/s]')
    if save == 1:
        plt.savefig(path_save_stat + 'time_track_norm.pdf')
    elif show == 1:
        plt.show()
    plt.clf()

def all_timetrack(x_values, *data, save_path, save=1, show=0):
    time_x = x_values
    path_save_stat = save_path
    omega_x_track, omega_x_err, omega_y_track, omega_y_err, omega_z_track, omega_z_err = data
    # plt.errorbar(time_x, omega_x_track, yerr=omega_x_err, errorevery=100, label='mean x-value of angular velocity', capsize=5, ecolor='#D4A368')
    # plt.errorbar(time_x, omega_y_track, yerr=omega_y_err, errorevery=100, label ='mean y-value of angular velocity', capsize=5, color='#D46899', ecolor='#57D49C')
    # plt.errorbar(time_x, omega_z_track, yerr=omega_z_err, errorevery=100, label ='mean z-value of angular velocity', capsize=5, color='#99D468', ecolor='#B579DE')
    plt.plot(time_x, omega_x_track,  label='mean x-value of angular velocity' )
    plt.plot(time_x, omega_y_track,  label ='mean y-value of angular velocity',  color='#D46899')
    plt.plot(time_x, omega_z_track,  label ='mean z-value of angular velocity',  color='#99D468')
    
    plt.legend()
    plt.title('Time profile of the angular velocity')
    plt.xlabel('time [days]')
    plt.ylabel('angular velocity [rad/s]')
    if save == 1:
        plt.savefig(path_save_stat + 'time_track_all.pdf')
    elif show == 1:
        plt.show()
    plt.clf()

def scatter_plot(x_values, *data, save_path, save=1, show=0):
    it_range = x_values
    path_save_scatter = save_path
    statistical_data, iteration, tx, ty, tz, w_x, w_y, w_z, omega_norm, rotation = data
    plt.scatter(it_range, tx, marker = '.', label = ' total tx data points')
    plt.scatter(it_range, ty, marker = 'x', label = 'total ty data points')
    plt.scatter(it_range, tz, marker = '.', label = 'total tz data points')
    plt.plot(it_range, [statistical_data['total tx'][0]] * iteration, label = ' total tx Mean value')
    plt.plot(it_range, [statistical_data['total ty'][0]] * iteration, label = 'total ty Mean value')
    plt.plot(it_range, [statistical_data['total tz'][0]] * iteration, label = 'total tz Mean value')
    plt.ylim([min(tz), max(tz)])
    plt.xlabel('experiment number')
    plt.ylabel(r'torque [Nm]')
    plt.title('Scatter plot of the final total torque components')
    plt.legend()
    if save == 1:
        plt.savefig(path_save_scatter + 'torque.pdf')
    elif show == 1:
        plt.show()
    plt.clf()

    plt.scatter( it_range,w_x, marker = '.', label = 'total wx data points')
    plt.scatter( it_range,w_y, marker = 'x', label = 'total wy data points')
    plt.scatter( it_range,w_z, marker = '.', label = 'total wz data points')
    plt.plot(it_range, [statistical_data['total wx'][0]] * iteration, label = ' wx Mean value')
    plt.plot(it_range, [statistical_data['total wy'][0]] * iteration, label = 'wy Mean value')
    plt.plot(it_range, [statistical_data['total wz'][0]] * iteration, label = 'wz Mean value')
    plt.ylim([min(w_z), max(w_z)])
    plt.xlabel('experiment')
    plt.ylabel(r'Angular velocity [$\frac{rad}{s}$]')
    plt.title('Scatter plot of the final total angular velocity components')
    plt.legend()
    if save == 1:
        plt.savefig(path_save_scatter + 'angular_velocity.pdf')
    elif show == 1:
        plt.show()
    plt.clf()
    
    

    plt.scatter(it_range, omega_norm, marker = 'x', label = ' angular speed (magnitude) data points')
    plt.plot(it_range, [statistical_data['|omega|'][0]] * iteration, label = ' angular speed (magnitude) Mean value')
    plt.ylim([min(omega_norm), max(omega_norm)])
    plt.xlabel('experiment number')
    plt.ylabel(r'Angular velocity [$frac{rad}{s}$]')
    plt.title('Scatter plot of the final angular velocity magnitude')
    plt.legend()
    if save == 1:
        plt.savefig(path_save_scatter + 'omega_norm.pdf')
    elif show == 1:
        plt.show()
    plt.clf()
    
    

    plt.scatter(it_range, rotation, marker = 'x', label = ' rotation /s (magnitude) data points')
    plt.plot(it_range, [statistical_data['rotation'][0]] * iteration, label = ' rotation /s (magnitude) Mean value')
    plt.ylim([min(rotation), max(rotation)])
    plt.xlabel('experiment number')
    plt.ylabel(r' full rotation [Hz]')
    plt.title('Scatter plot of the final rotation of the CubeSat')
    plt.legend()
    if save == 1:
        plt.savefig(path_save_scatter + 'rotation.pdf')
    elif show == 1:
        plt.show()
    plt.clf()

def hist_plot(nbins, *data, save_path, save=1, show=0):
    path_save_stat = save_path
    w_x, w_y, w_z, omega_norm, w_x_stat, w_y_stat, w_z_stat, omega_stat = data
    bin_centres_z, hist_fit_z, start_z, end_z = hist_ply(w_z, w_z_stat[0], w_z_stat[1], nbins)
    bin_centres_y, hist_fit_y, start_y, end_y = hist_ply(w_y, w_y_stat[0], w_y_stat[1], nbins)
    bin_centres_x, hist_fit_x, start_x, end_x = hist_ply(w_x, w_x_stat[0], w_x_stat[1], nbins)
    bin_centres, hist_fit, start, end = hist_ply(omega_norm, omega_stat[0], omega_stat[1], nbins)

    plt.hist(w_x, bins=nbins, label='x-components of angular velocity', color='#578abd')
    plt.plot(bin_centres_x, hist_fit_x, label='Gaussian fit', color='r')
    plt.plot([w_x_stat[0]] * int(max(hist_fit_x)), range(int(max(hist_fit_x))), label='mean value', color='#22d800') 
    plt.plot([w_x_stat[0] + w_x_stat[1]] * int(max(hist_fit_x)), range(int(max(hist_fit_x))), label='standard deviation', color ='#FF8C00')
    plt.plot([w_x_stat[0] - w_x_stat[1]] * int(max(hist_fit_x)), range(int(max(hist_fit_x))), color='#FF8C00')
    plt.legend()
    plt.title('x-component of angular velocity distribution')
    plt.xlabel('Angular velocity bins [rad/s]')
    plt.ylabel('Counts per bin [-]')
    if save == 1:
        plt.savefig(path_save_stat + 'angular_x_stat.pdf')
    elif show == 1:
        plt.show()
    plt.clf()
    
    
    plt.hist(w_y, bins=nbins, label='y-components of angular velocity', color='#578abd')
    plt.plot(bin_centres_y, hist_fit_y, label='Gaussian fit', color='r')
    plt.plot([w_y_stat[0]] * int(max(hist_fit_y)), range(int(max(hist_fit_y))), label='mean value', color='#22d800') 
    plt.plot([w_y_stat[0] + w_y_stat[1]] * int(max(hist_fit_y)), range(int(max(hist_fit_y))), label='standard deviation', color ='#FF8C00')
    plt.plot([w_y_stat[0] - w_y_stat[1]] * int(max(hist_fit_y)), range(int(max(hist_fit_y))), color='#FF8C00')
    plt.legend()
    plt.title('y-component of angular velocity distribution')
    plt.xlabel('Angular velocity bins [rad/s]')
    plt.ylabel('Counts per bin [-]')
    if save == 1:
        plt.savefig(path_save_stat + 'angular_y_stat.pdf')
    elif show == 1:
        plt.show()
    plt.clf() 

    
    

    plt.hist(w_z, bins=nbins, label='z-components of angular velocity', color='#578abd')
    plt.plot(bin_centres_z, hist_fit_z, label='Gaussian fit', color='r')
    plt.plot([w_z_stat[0]] * int(max(hist_fit_z)), range(int(max(hist_fit_z))), label='mean value', color='#22d800') 
    plt.plot([w_z_stat[0] + w_z_stat[1]] * int(max(hist_fit_z)), range(int(max(hist_fit_z))), label='standard deviation', color ='#FF8C00')
    plt.plot([w_z_stat[0] - w_z_stat[1]] * int(max(hist_fit_z)), range(int(max(hist_fit_z))), color='#FF8C00')
    plt.legend()
    plt.title('z-component of angular velocity distribution')
    plt.xlabel('Angular velocity bins [rad/s]')
    plt.ylabel('Counts per bin [-]')
    if save == 1:
        plt.savefig(path_save_stat + 'angular_z_stat.pdf')
    elif show == 1:
        plt.show()
    plt.clf()
    

    plt.hist(omega_norm, bins=nbins, label='norm of angular velocity', color='#578abd')
    plt.plot(bin_centres, hist_fit, label='Gaussian fit', color='r')
    plt.plot([omega_stat[0]] * int(max(hist_fit)), range(int(max(hist_fit))), label='mean value', color='#22d800') 
    plt.plot([omega_stat[0] + omega_stat[1]] * int(max(hist_fit)), range(int(max(hist_fit))), label='standard deviation', color ='#FF8C00')
    plt.plot([omega_stat[0] - omega_stat[1]] * int(max(hist_fit)), range(int(max(hist_fit))), color='#FF8C00')
    plt.legend()
    plt.title('Norm of angular velocity distribution')
    plt.xlabel('Angular velocity bins [rad/s]')
    plt.ylabel('Counts per bin [-]')
    if save == 1:
        plt.savefig(path_save_stat + 'angular_norm_stat.pdf')
    elif show == 1:
        plt.show()
    plt.clf()

def plot_histogram(nbins, data, stat_data, title_str, x_label, name, save_path, save=1, show=0):

    bin_centres, hist_fit, start, end = hist_ply(data, stat_data[0], stat_data[1], nbins)
    
    #fit the gaussian curve
    f_x = []
    for i in bin_centres:
        f_x.append(gauss(i, max(hist_fit), stat_data[0], stat_data[1]))

    plt.hist(data, bins=nbins, label='real count', color='#578abd')
    plt.plot(bin_centres, f_x, label='Gaussian fit', color='r')
    plt.plot([stat_data[0]] * int(max(hist_fit)), range(int(max(hist_fit))), label='mean value', color='#22d800') 
    plt.plot([stat_data[0] + stat_data[1]] * int(max(hist_fit)), range(int(max(hist_fit))), label='standard deviation', color ='#FF8C00')
    plt.plot([stat_data[0] - stat_data[1]] * int(max(hist_fit)), range(int(max(hist_fit))), color='#FF8C00')
    plt.legend()
    plt.title(title_str)
    plt.xlabel(x_label)
    plt.ylabel('Counts per bin [-]')
    if save == 1:
        plt.savefig(save_path + '{}.pdf'.format(name))
    elif show == 1:
        plt.show()
    plt.clf()
def improvement(x_data,data_baseline, data, path_save_stat, show=0, save=1):
    improv = []
    for i in range(len(x_data)):
        improv = np.append(improv, ((data_baseline[i] - data[i]) / data_baseline[i])*100)
    
    plt.plot(x_data, improv, label='improvement percentage')
    plt.title(r'Improvement of Angular Velocity Norm')
    plt.xlabel('Time [days]')
    plt.ylabel('Percentage [%]')
    plt.legend()
    if show == 1:
        plt.show()
    elif save == 1:
        plt.savefig(path_save_stat)


def plot_single_timetrack(time, data, savepath, title_str, label_str, show=0, save=1):
    plt.plot(time, data, color='#FFA500', label=label_str) #D4A368
    plt.title(title_str)
    plt.xlabel('Time [days]')
    plt.ylabel('Angular velocity [rad/s]')
    plt.legend()
    if save == 1:
        plt.savefig(savepath)
    if show ==1:
        plt.show() 
    plt.clf()


def plot_grid_values(value_dict, grid, cbar_unit, title_str, savepath, quadrant_sep=True, save=1, show=0):
    x = []
    y = []
    tb = []
    markersize = []
    n = 0
    for i in value_dict.keys():
        x.append(grid.pixelPosition[i][1])
        y.append(grid.pixelPosition[i][2])
        tb.append(100*value_dict[i]/(86400*365))
        markersize.append(tb[n] * 20)
        n+=1


    plt.scatter(x, y, c=tb) #colour; c=tb

    if quadrant_sep == True:
        plt.plot([0]*len(y), y, c='red')
        plt.plot(x, [0]*len(x), c='red')
    cbar = plt.colorbar()
    cbar.set_label(cbar_unit)
    plt.xlabel('y-position [m]')
    plt.ylabel('z-position [m]')
    plt.title(title_str)
    if save == 1:
        plt.savefig(savepath)
    if show ==1:
        plt.show() 
    plt.clf()


def centerQuadrantPlot(value_array, dict_list, cbar_unit, title_str, savepath, quadrant_sep=True, save=1, show=0):
    x = []
    y = []

    x_i = []
    y_i = []
    
    for dict in dict_list:
        x_i = []
        y_i = []
        for pix in dict.keys():

            x_i.append(dict[pix][:][1])
            y_i.append(dict[pix][:][2])

        x.append(np.mean(x_i))
        y.append(np.mean(y_i))

    scale = 100*value_array/(86400*365)
    plt.scatter(x, y, s=2500*scale, c=scale) #colour; c=tb


    if quadrant_sep == True:
        plt.plot([0]*100, np.linspace(-0.03, 0.03, 100), c='red')
        plt.plot(np.linspace(-0.03, 0.03, 100), [0]*100, c='red')
    cbar = plt.colorbar()
    cbar.set_label(cbar_unit)
    plt.xlabel('y-position [m]')
    plt.ylabel('z-position [m]')

    plt.ylim((-0.032, 0.032))
    plt.xlim((-0.032, 0.032))
    plt.title(title_str)
    if save == 1:
        plt.savefig(savepath)
    if show ==1:
        plt.show() 
    plt.clf()
