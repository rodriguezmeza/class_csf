#
# Execute:
# $ python getdist_plots_no-nu_csf_no-flat.py > getdist_no-nu_csf_no-flat.txt
#

# Export the results to GetDist
from cobaya import load_samples
chaindir1 = "chains_no-flat/chain"
legend1 = "BAO-SN-Planck"
color1 = 'orange'
ls1 = '-'
line1 = {'ls':ls1, 'lw':1.0, 'color':color1}
folder = "results_bao_dr2_no-flat/"
huevos_triangle = folder + "huevos_csf_bao_dr2.pdf"
#
huevos_Om_OL_2d = folder + "huevos_csf_bao_dr2_Om_OL_2d.pdf"
huevos_Om_OK_2d = folder + "huevos_csf_bao_dr2_Om_OK_2d.pdf"
huevos_H0_Om_2d = folder + "huevos_csf_bao_dr2_H0_Om_2d.pdf"
huevos_Om_1d = folder + "huevos_emt_csf_dr2_Om_1d.pdf"
covariance_file_1 = folder + "covmat_csf_bao_dr2.txt"

parameters = ["Omega_m",  "Omega_Lambda", "H0_CSF", "omega_b",
            "Omega_k", "QCSF"]
parameters_Om_OL = ["Omega_m",  "Omega_Lambda"]
parameters_Om_OK = ["Omega_m",  "Omega_k"]
parameters_H0_Om = ["H0_CSF",  "Omega_m"]

legends = [legend1]
lineargs=[line1]
lineargs_triangle=[line1]
#legendloc='upper right'
legendloc=''
contourcolors=[color1]

#DESI-VI 2024 (2404.03002):
# LCDM+Omega_k : DESI + BBN + theta_* :
# Omega_m = 0.284 +- 0.020
# Omega_m (DESI-CMB) = 0.3069 +- 0.0050
# Omega_Lambda = 0.651 (+ 0.068 -0.057)
# H0 = 68.52 +- 0.69
# H0 (DESI-CMB, flat) = 67.97 +- 0.38
# Omega_k (10^3) = 0.3 +4.8 -5.4
# Omega_k (CMB) = -0.0102 +-0.0054
# Omega_k (DESI-CMB) = 0.00224 +-0.0016
#
# LCDM+Omega_k : DESI + BBN :
# Omega_m (only DESI) = 0.284 +- 0.020
# Omega_m (DESI-CMB) = 0.3049 +- 0.0051
# H0 (DESI-CMB) = 68.51 +- 0.52
# Omega_k (10^3) = 65 +68 -78
# Omega_k (10^3) (DESI-CMB) = 2.4 +-1.6 -> 2.4x10^-3 -> 0.0024
# Omega_Lambda (only DESI) = 0.651 (+ 0.068 -0.057)
# Omega_Lambda (DESI-CMB, 1-Om-OK) = 0.6927
#
# sigma8 (68%, CMB, lensing + BAO) = 0.811 +-0.019
# sigma8 (plik,CamSpec, combined) = 0.8101 +-0.0061
# sigma8 (%68, TT,TE,EE+lowE+lensing+BAO, table 2 planck-paper) = 0.8102 +-0.0060
#
redMarkersLabels = {'Omega_m': 0.3142, 'Omega_Lambda':0.6858, 'H0':67.39,
    'omega_b': 0.02229, 'Omega_k': 0.00, 'sigma8':0.8091, 'nCSF':0}
marker_args={'ls':'--','lw':0.75, 'color':'gray'}

# This is cobaya load samples
#gd_sample1 = load_samples(chaindir1, to_getdist=True)

from getdist import loadMCSamples, chains, MCSamples
#
# 'mult_bias_correction_order':1,0
# 'smooth_scale_2D':0.55,0.3
# 'smooth_scale_1D':0.55,0.3
#
gd_sample1 = loadMCSamples(chaindir1, settings={'ignore_rows':0.3,
    'mult_bias_correction_order':1, 'smooth_scale_2D':0.55, 'smooth_scale_1D':0.55})

pranges = gd_sample1.setRanges({'Omega_m':[0.,0.5]})
pranges = gd_sample1.setRanges({'Omega_Lambda':[0.3,0.8]})
pranges = gd_sample1.setRanges({'QCSF':[-0.75,0.75]})

# %matplotlib inline  # uncomment if running from the Jupyter notebook
import getdist.plots as gdplt
import matplotlib.pyplot as plt
import numpy as np

legend_fontsize = 9
axes_fontsize = 10
lab_fontsize = 15

samples_bao_smnu=[gd_sample1]
legends_triangle_bao_smnu = [legend1]
contourcolors_triangle_bao_smnu=[color1]
parameters_triangle = ["Omega_m", "H0_CSF", "Omega_Lambda",
            "Omega_k", "QCSF", "T_cmb_CSF"]
gdplot = gdplt.get_subplot_plotter()
gdplot.settings.legend_fontsize = 20
gdplot.settings.axes_fontsize = 20
gdplot.settings.lab_fontsize = 30
gdplot.triangle_plot(samples_bao_smnu,
    parameters_triangle, filled=True,
    legend_labels=legends_triangle_bao_smnu, legend_loc=legendloc,
    contour_colors=contourcolors_triangle_bao_smnu, line_args=lineargs_triangle,
    markers=redMarkersLabels,
    marker_args=marker_args,
    fontsize=20,
    )
gdplot.export(huevos_triangle)


g=gdplt.get_single_plotter(width_inch=4, ratio=0.9)
g.settings.legend_fontsize = legend_fontsize
g.settings.axes_fontsize = axes_fontsize
g.settings.lab_fontsize = lab_fontsize
g.settings.figure_legend_frame = False
g.plot_2d([gd_sample1], parameters_Om_OL, filled=True,
    colors=[color1],
    lims=[0.25,0.4,0.55,0.85]
    )
g.add_legend([legend1], legend_loc='lower right')
x=np.linspace(0, 1, 100)
plt.plot(x,1-x, color='gray',ls='--', lw=0.75)
x=np.linspace(0, 1, 100)
plt.plot(x,0.6858+0*x, color='gray',ls='--', lw=0.75)
y=np.linspace(0, 1, 100)
plt.plot(0.3142+0*y,y, color='gray',ls='--', lw=0.75)
g.export(huevos_Om_OL_2d)

g=gdplt.get_single_plotter(width_inch=4, ratio=0.9)
g.settings.legend_fontsize = legend_fontsize
g.settings.axes_fontsize = axes_fontsize
g.settings.lab_fontsize = lab_fontsize
g.settings.figure_legend_frame = False
g.plot_2d([gd_sample1], parameters_Om_OK, filled=True,
    colors=[color1],
    lims=[0.2,0.4,-0.2,0.2], markers={'Omega_k':0.003}
    )
g.add_legend([legend1], legend_loc='upper right')
x=np.linspace(0, 1, 100)
plt.plot(x,0*x, color='gray',ls='--', lw=0.75)
y=np.linspace(-0.5, 0.5, 100)
plt.plot(0.3142+0*y,y, color='gray',ls='--', lw=0.75)
g.export(huevos_Om_OK_2d)


#B Om vs H0 2D
g=gdplt.get_single_plotter(width_inch=4, ratio=0.9)
g.settings.legend_fontsize = legend_fontsize
g.settings.axes_fontsize = axes_fontsize
g.settings.lab_fontsize = lab_fontsize
g.settings.figure_legend_frame = False
g.plot_2d([gd_sample1], parameters_H0_Om, filled=True,
    colors=[color1],
    lims=[60,80,0.2,0.4]
    )
g.add_legend([legend1], legend_loc='upper right')
x=np.linspace(0, 100, 100)
plt.plot(x,0.3142+0*x, color='gray',ls='--', lw=0.75)
y=np.linspace(0, 1, 100)
plt.plot(67.39+0*y,y, color='gray',ls='--', lw=0.75)
g.export(huevos_H0_Om_2d)
#E

g=gdplt.get_single_plotter(width_inch=4, ratio=0.9)
g.settings.legend_fontsize = legend_fontsize
g.settings.axes_fontsize = axes_fontsize
g.settings.lab_fontsize = lab_fontsize
g.settings.figure_legend_frame = False
g.plot_1d([gd_sample1],
    "Omega_m", colors=[color1],
        lims=[0.225,0.4],
        normalized=True
        )
y=np.linspace(0, 100, 100)
plt.plot(0.3142+0*y,y, color='gray',ls='--', lw=0.75)
g.add_legend([legend1], legend_loc='upper right')
g.export(huevos_Om_1d)

# Many statistic analyses

#import sys, io
#stdout = sys.stdout
#sys.stdout = io.StringIO()
#
#%%capture cap --no-stderr

# Analyzing
mean = gd_sample1.getMeans()[:6]
covmat = gd_sample1.getCovMat().matrix[:6, :6]
print("Mean:")
print(mean)
print("Covariance matrix:")
print(covmat)
print()


print()
print("Gelman-Rubin analysis::")
print()
print("bao_dr2")
print(gd_sample1.getGelmanRubin())
print(gd_sample1.getGelmanRubinEigenvalues())
print()

print()
print("PCA analysis::")
print()
print(gd_sample1.PCA(parameters))
print()
#

print()
print("Covariance analysis::")
pcovmat = gd_sample1.getCovMat()
print()
print('Covariance matrix: ',pcovmat.matrix)
print()
pcovmat.saveToFile(covariance_file_1)
print()

# Many other things you can do besides plot, e.g. get latex
# Default limits are 1: 68%, 2: 95%, 3: 99% probability enclosed
# See  https://getdist.readthedocs.io/en/latest/analysis_settings.html
# and examples for below for changing analysis settings
# (e.g. 2hidh limits, and how they are defined)

print()
print("getInlineLatex analysis::")
print()
print(gd_sample1.getInlineLatex('Omega_m',limit=2))
print(gd_sample1.getInlineLatex('Omega_Lambda',limit=2))
print(gd_sample1.getInlineLatex('H0',limit=2))
print(gd_sample1.getInlineLatex('omega_b',limit=2))
print(gd_sample1.getInlineLatex('Omega_k',limit=2))
print(gd_sample1.getInlineLatex('H0_CSF',limit=2))
print(gd_sample1.getInlineLatex('T_cmb_CSF',limit=2))
print(gd_sample1.getInlineLatex('QCSF',limit=2))
print()

print()
print("getTable analysis::")
print()
print(gd_sample1.getTable().tableTex())
print()


# results from multiple chains
parameters_ResultsTable = ["logA", "n_s", "Omega_k", "H0_CSF", "omega_b", "omega_cdm", "tau_reio", "Omega_m",  "Omega_Lambda", "H0", "omega_b",
            "QCSF",
            "A_s"]
print()
print()
print("Results table::")
print()
results = [gd_sample1]
title1 = 'BAO DR2'
titles = [title1]
from getdist.types import ResultTable
print(ResultTable(ncol=1,results=results,
                 paramList=parameters_ResultsTable, limit=1, titles=titles).tableTex())
print()

# if samples have likelihood values, can also get best fit sample and extremal values of N-D confidence region
# Note in high dimensions best-fit sample is likely a long way from the best fit; N-D limits also often MC-noisy
print()
print("Stats table::")
print()
print(gd_sample1.getLikeStats())
print()


print()
print("confidence extrema table::")
print()

print('Omega_m 95% n-D confidence extrema:',
    gd_sample1.paramNames.parWithName('Omega_m').ND_limit_bot[1],
    gd_sample1.paramNames.parWithName('Omega_m').ND_limit_top[1])
print()

print('Omega_Lambda 95% n-D confidence extrema:',
    gd_sample1.paramNames.parWithName('Omega_Lambda').ND_limit_bot[1],
    gd_sample1.paramNames.parWithName('Omega_Lambda').ND_limit_top[1])
print()

print('H0_CSF 95% n-D confidence extrema:',
    gd_sample1.paramNames.parWithName('H0_CSF').ND_limit_bot[1],
    gd_sample1.paramNames.parWithName('H0_CSF').ND_limit_top[1])
print()

print('T_cmb_CSF 95% n-D confidence extrema:',
    gd_sample1.paramNames.parWithName('T_cmb_CSF').ND_limit_bot[1],
    gd_sample1.paramNames.parWithName('T_cmb_CSF').ND_limit_top[1])
print()

print('QCSF 95% n-D confidence extrema:',
    gd_sample1.paramNames.parWithName('QCSF').ND_limit_bot[1],
    gd_sample1.paramNames.parWithName('QCSF').ND_limit_top[1])
print()


#with open(getdist_all_info_txt_file, 'w') as f:
#    f.write(cap.stdout)
