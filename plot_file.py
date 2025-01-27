import scipy.constants as sc
import scipy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from os import path

#convert separation in AU to period in years
def period(total_mass,separation):
    return(separation ** 3/total_mass) ** (1/2)

AU = sc.au

#sndata = np.genfromtxt('../LISA_sensitivities410.txt')
sndata = np.genfromtxt('../noise_characterictic_strain.txt')
fn = sndata[:,0]
sn4 = sndata[:,1]
#sn10 = sndata[:,5]

#to interpolate noise curves at any frequency between max and min in the file use:
noise_int4 = scipy.interpolate.interp1d(fn,sn4)
#noise_int10 = scipy.interpolate.interp1d(fn,sn10)

T4 = 4*3.15e7 # 4yr in seconds representing LISA mission nominal lifetime 
#T10 = 10*3.15e7 #10yr extended LISA mission

triple = np.load('../triple_starp_data.npy')
triple_starp_data = pd.DataFrame(triple, columns=['m1s', 'm2s', 'a', 'age_starparticle', 'mass_starparticle',  'distance_starparticle', 'system_index','gal_index','st1','st2','e','x','y','z','sys_time'])
print("triple read",len(triple))

binary = np.load('../binary_HR_starp_data.npy')
binary_starp_data = pd.DataFrame(binary, columns=['m1s', 'm2s', 'a', 'age_starparticle', 'mass_starparticle',  'distance_starparticle', 'system_index','gal_index','st1','st2','e','x','y','z','sys_time'])
print("binary read",len(binary))

cir_triple_starp_data = triple_starp_data[triple_starp_data.e < 0.5]
ecc_triple_starp_data = triple_starp_data[triple_starp_data.e > 0.5]

m1_cir = cir_triple_starp_data.m1s
m2_cir = cir_triple_starp_data.m2s
total_mass_cir = m1_cir+ m2_cir
a_cir = cir_triple_starp_data.a
e_cir = cir_triple_starp_data.e
index_cir = cir_triple_starp_data['system_index']
st1_cir = cir_triple_starp_data['st1']
st2_cir = cir_triple_starp_data['st2']

r_cir = cir_triple_starp_data.distance_starparticle*1000*3.086e+16 
age_cir = cir_triple_starp_data.age_starparticle
gal_index_cir = cir_triple_starp_data['gal_index']

gw_f_cir = 2/(period(total_mass_cir,a_cir)*365*24*60*60)
mchirp_cir = (m1_cir*m2_cir)**0.6 / (m1_cir+m2_cir)**0.2
Amp_gw_cir = (2*(sc.G*mchirp_cir*2e30)**(5./3.)*(gw_f_cir*np.pi)**(2./3.))/(r_cir*sc.c**4)

m1_ecc = ecc_triple_starp_data.m1s
m2_ecc = ecc_triple_starp_data.m2s
total_mass_ecc = m1_ecc+m2_ecc
a_ecc = ecc_triple_starp_data.a
e_ecc = ecc_triple_starp_data.e
index_ecc = ecc_triple_starp_data['system_index']
st1_ecc = ecc_triple_starp_data['st1']
st2_ecc = ecc_triple_starp_data['st2']

r_ecc = ecc_triple_starp_data.distance_starparticle*1000*3.086e+16 
age_ecc = ecc_triple_starp_data.age_starparticle
gal_index_ecc = ecc_triple_starp_data['gal_index']



f_orb = 1/(period(total_mass_ecc,a_ecc)*365*24*60*60)
f_bur = 2 * f_orb * ((1-e_ecc) ** (-3/2))
mchirp_ecc = (m1_ecc*m2_ecc)**0.6 / (m1_ecc+m2_ecc)**0.2
#Amp_bur = np.sqrt(32/5) * m1_ecc*m2_ecc / (r_ecc*(1-e_ecc)) 
eta  =  4*m1_ecc*m2_ecc / ((total_mass_ecc)**2)
Amp_bur = 7.65e-21 * eta * ((total_mass_ecc/20)**(5/3)) * ((f_bur/3.16e-3)**(2/3)) * ((r_ecc/(8*1000*3.086e+16))**-1)
t_p = (1-e_ecc)**(3/2) * (period(total_mass_ecc,a_ecc)*365*24*60*60)

m1 = triple_starp_data.m1s
m2 = triple_starp_data.m2s
total_mass = m1+m2
r = triple_starp_data.distance_starparticle*1000*3.086e+16 
a = triple_starp_data.a
e = triple_starp_data.e
age = triple_starp_data.age_starparticle
index = triple_starp_data['system_index']
gal_index = triple_starp_data['gal_index']
st1 = triple_starp_data['st1']
st2 = triple_starp_data['st2']
mchirp = (m1*m2)**0.6 / (m1+m2)**0.2

gw_f = pd.concat([f_bur, gw_f_cir], ignore_index=True)
Amp_gw = pd.concat([Amp_bur, Amp_gw_cir], ignore_index=True)

cir_binary_starp_data = binary_starp_data[binary_starp_data.e < 0.5]
ecc_binary_starp_data = binary_starp_data[binary_starp_data.e > 0.5]

m1_cir_binary = cir_binary_starp_data.m1s
m2_cir_binary = cir_binary_starp_data.m2s
total_mass_cir_binary = m1_cir_binary+m2_cir_binary
a_cir_binary = cir_binary_starp_data.a
e_cir_binary = cir_binary_starp_data.e
index_cir_binary = cir_binary_starp_data['system_index']
st1_cir_binary = cir_binary_starp_data['st1']
st2_cir_binary = cir_binary_starp_data['st2']

r_cir_binary = cir_binary_starp_data.distance_starparticle*1000*3.086e+16 
age_cir_binary = cir_binary_starp_data.age_starparticle
gal_index_cir_binary = cir_binary_starp_data['gal_index']

gw_f_cir_binary = 2/(period(total_mass_cir_binary,a_cir_binary)*365*24*60*60)
mchirp_cir_binary = (m1_cir_binary*m2_cir_binary)**0.6 / (m1_cir_binary+m2_cir_binary)**0.2
Amp_gw_cir_binary = (2*(sc.G*mchirp_cir_binary*2e30)**(5./3.)*(gw_f_cir_binary*np.pi)**(2./3.))/(r_cir_binary*sc.c**4)

m1_ecc_binary = ecc_binary_starp_data.m1s
m2_ecc_binary = ecc_binary_starp_data.m2s
total_mass_ecc_binary = m1_ecc_binary+m2_ecc_binary
a_ecc_binary = ecc_binary_starp_data.a
e_ecc_binary = ecc_binary_starp_data.e
index_ecc_binary = ecc_binary_starp_data['system_index']
st1_ecc_binary = ecc_binary_starp_data['st1']
st2_ecc_binary = ecc_binary_starp_data['st2']

r_ecc_binary = ecc_binary_starp_data.distance_starparticle*1000*3.086e+16 
age_ecc_binary = ecc_binary_starp_data.age_starparticle
gal_index_ecc_binary = ecc_binary_starp_data['gal_index']



f_orb_binary = 1/(period(total_mass_ecc_binary,a_ecc_binary)*365*24*60*60)
f_bur_binary = 2 * f_orb_binary * ((1-e_ecc_binary) ** (-3/2))
mchirp_ecc_binary = (m1_ecc_binary*m2_ecc_binary)**0.6 / (m1_ecc_binary+m2_ecc_binary)**0.2
#Amp_bur_binary = np.sqrt(32/5) * m1_ecc_binary*m2_ecc_binary / (r_ecc_binary*(1-e_ecc_binary))
eta_binary  =  4*m1_ecc_binary*m2_ecc_binary / ((total_mass_ecc_binary)**2)
Amp_bur_binary = 7.65e-21 * eta_binary * ((total_mass_ecc_binary/20)**(5/3)) * ((f_bur_binary/3.16e-3)**(2/3)) * ((r_ecc_binary/(8*1000*3.086e+16))**-1)
t_p_binary = (1-e_ecc_binary)**(3/2) * (period(total_mass_ecc_binary,a_ecc_binary)*365*24*60*60)

m1_binary = binary_starp_data.m1s
m2_binary = binary_starp_data.m2s
total_mass_binary = m1_binary+m2_binary
r_binary = binary_starp_data.distance_starparticle*1000*3.086e+16 
a_binary = binary_starp_data.a
e_binary = binary_starp_data.e
age_binary = binary_starp_data.age_starparticle
index_binary = binary_starp_data['system_index']
gal_index_binary = binary_starp_data['gal_index']
st1_binary = binary_starp_data['st1']
st2_binary = binary_starp_data['st2']
mchirp_binary = (m1_binary*m2_binary)**0.6 / (m1_binary+m2_binary)**0.2

gw_f_binary = pd.concat([f_bur_binary, gw_f_cir_binary], ignore_index=True)
Amp_binary = pd.concat([Amp_bur_binary, Amp_gw_cir_binary], ignore_index=True)

snr4 = 0.5*Amp_gw*np.sqrt(T4)/noise_int4(gw_f)
snr4_binary = 0.5*Amp_binary*np.sqrt(T4)/noise_int4(gw_f_binary)

# Filter values where SNR > 7
filtered_indices = snr4 > 7
filtered_indices_binary = snr4_binary > 7

filtered_Amp_gw = Amp_gw[filtered_indices]
filtered_gw_f = gw_f[filtered_indices]
filtered_m1 = m1[filtered_indices]
filtered_m2 = m2[filtered_indices]
filtered_mchirp = mchirp[filtered_indices]
filtered_e = e[filtered_indices]

filtered_m1_binary = m1_binary[filtered_indices_binary]
filtered_m2_binary = m2_binary[filtered_indices_binary]
filtered_e_binary = e_binary[filtered_indices_binary]
filtered_mchirp_binary = mchirp_binary[filtered_indices]
filtered_Amp_gw_binary = Amp_binary[filtered_indices_binary]
filtered_gw_f_binary = gw_f_binary[filtered_indices_binary]

# Combine the two Series into a DataFrame
triple_st = pd.DataFrame({'st1': st1, 'st2': st2})

# Ensure each combination is sorted, so (10, 11) and (11, 10) are treated the same
triple_st_sorted = triple_st.apply(lambda x: sorted([x['st1'], x['st2']]), axis=1)
triple_st_sorted = pd.DataFrame(triple_st_sorted.tolist(), columns=['st1', 'st2'])

# Count occurrences of each unique combination
triple_combination_counts = triple_st_sorted.groupby(['st1', 'st2']).size().reset_index(name='count')

# Display the result
print(triple_combination_counts)

# Combine the two Series into a DataFrame
binary_st = pd.DataFrame({'st1_binary': st1_binary, 'st2_binary': st2_binary})

# Ensure each combination is sorted, so (10, 11) and (11, 10) are treated the same
binary_st_sorted = binary_st.apply(lambda x: sorted([x['st1_binary'], x['st2_binary']]), axis=1)
binary_st_sorted = pd.DataFrame(binary_st_sorted.tolist(), columns=['st1_binary', 'st2_binary'])

# Count occurrences of each unique combination
binary_combination_counts = binary_st_sorted.groupby(['st1_binary', 'st2_binary']).size().reset_index(name='count')

# Display the result
print(binary_combination_counts)

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc

#Enable LaTeX for all text in the plots
rc('text', usetex=True)
rc('font', family='serif')

# Customizable parameters
fontsize = 28
bar_width = 0.35
figsize = (12, 7)
title = 'Comparison of Combination Counts'
title_fontsize = 16
title_weight = 'bold'
xlabel = 'Labels'
xlabel_fontsize = fontsize
ylabel = 'Counts'
ylabel_fontsize = fontsize
xticks_fontsize = fontsize
yticks_fontsize = fontsize
tick_direction = 'inout'
tick_length = 6
tick_width = 2
grid_linestyle = '--'
grid_alpha = 0.7
grid_linewidth = 0.8
bar_edgecolor = 'black'
bar_linewidth = 1.5
bar_color1 = 'green'
bar_alpha1 = 0.6
label1 = 'Triple'
bar_color2 = 'purple'
bar_alpha2 = 0.6
label2 = 'Binary'
legend_fontsize = 0.9*fontsize
legend_loc = 'upper right'
figure_edgecolor = 'black'  # Figure border color
figure_linewidth = 6  # Figure border thickness
axes_border_thickness = 2

# Data
#labels = ['He-He', 'He-CO', 'He-ONe', 'CO-CO', 'CO-ONe','ONe-ONe']
#values_triple = triple_combination_counts['count']
#values_binary = binary_combination_counts['count']
labels = ['He-He', 'He-CO', 'He-ONe', 'CO-CO', 'CO-ONe', 'ONe-ONe']

values_triple = [2632919, 3866966, 314439, 345812, 101990, 130]
values_binary = [1668050, 1901346, 85279, 194019, 59700] + [0]  # Append 0 to match the length


# Set up bar positions
x = range(len(labels))

# Create the bar chart
fig, ax = plt.subplots(figsize=figsize)

# Plot bars with customizable border thickness
bars1 = ax.bar([p - bar_width/2 for p in x], values_triple, width=bar_width, color=bar_color1, alpha=bar_alpha1, 
               label=label1, edgecolor=bar_edgecolor, linewidth=bar_linewidth)
bars2 = ax.bar([p + bar_width/2 for p in x], values_binary, width=bar_width, color=bar_color2, alpha=bar_alpha2, 
               label=label2, edgecolor=bar_edgecolor, linewidth=bar_linewidth)

ax.set_yscale('log')

# Add x-axis tick labels
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=xticks_fontsize)

#yticklabels  = [

for label in ax.get_yticklabels():
    label.set_fontsize(yticks_fontsize)

ax.tick_params(axis='x', direction=tick_direction, length=tick_length, width=tick_width, top=True, bottom=True)
ax.tick_params(axis='y', direction=tick_direction, length=tick_length, width=tick_width, left=True, right=True)


# Show grid lines for better readability
ax.yaxis.grid(True, linestyle=grid_linestyle, alpha=grid_alpha, linewidth=grid_linewidth)

# Add value labels on top of bars
for bar in bars1:
    height = bar.get_height()

for bar in bars2:
    height = bar.get_height()

# Add legend
ax.legend(fontsize=legend_fontsize, loc=legend_loc)

# Adjust the border thickness of the axes
for spine in ax.spines.values():
    spine.set_linewidth(axes_border_thickness)

# Adjust layout to make room for rotated labels
plt.tight_layout()
plt.savefig('st1_binary_triple.pdf')

plt.show()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

#Enable LaTeX for all text in the plots
rc('text', usetex=True)
rc('font', family='serif')

# Customizable parameters
background_color = 'white'  # Light Gray 2 from Google Sheets
border_thickness = 2
line_thickness = 2
fontsize = 28
ticksize = 24
# Create the plot
fig, ax = plt.subplots(figsize=(10, 8))
fig.patch.set_facecolor(background_color)
ax.set_facecolor(background_color)

# Plot scatter plots
marker_size = 150
# Plot scatter points
ax.scatter(gw_f, Amp_gw * np.sqrt(T4), color='grey', marker='*', s=marker_size, label='Unresolved DWD (Triple )',rasterized=True)
ax.scatter(f_bur, Amp_bur * np.sqrt(t_p), color='black', marker='*', s=marker_size, label='Eccentric (Triple)',rasterized=True)
ax.scatter(filtered_gw_f, filtered_Amp_gw * np.sqrt(T4), color='green', marker='*', s=marker_size, label='Resolved DWD (Triple)',rasterized=True)
ax.scatter(filtered_gw_f_binary, filtered_Amp_gw_binary * np.sqrt(T4), color='purple', label='Resolved DWD (Binary)',rasterized=True)

# Plot LISA noise curves
ax.loglog(fn, sn4, 'k-', linewidth=line_thickness, label='LISA noise')

# Set plot parameters
ax.tick_params(axis='both', which='major', labelsize=ticksize, direction='inout', length=6)
ax.set_ylabel('$\mathrm{h_{c}}$', fontsize=fontsize)  # y-axis label
ax.set_xlabel('$\mathrm{f_{gw}(Hz)}$', fontsize=fontsize)  # x-axis label
ax.set_xlim(1e-4, 1e-1)  # x-axis limits
ax.legend(ncol=3, fontsize=18, bbox_to_anchor=(1.125, 1.25))

# Customize plot border
for spine in ax.spines.values():
    spine.set_linewidth(border_thickness)

# Set ticks on all sides
ax.tick_params(top=True, bottom=True, left=True, right=True)

# Adjust margins to prevent labels from being cut off
plt.subplots_adjust(left=0.18, right=0.915, top=0.7, bottom=0.18)

# Save the plot
plt.savefig('triple_hf1.png')
plt.savefig('triple_hf1.pdf')
plt.show()
print("triple plot done")

primary_mass = np.maximum(triple_starp_data.m1s, triple_starp_data.m2s)
filtered_primary_mass = np.maximum(filtered_m1, filtered_m2)

primary_mass_binary = np.maximum(binary_starp_data.m1s, binary_starp_data.m2s)
filtered_primary_mass_binary = np.maximum(filtered_m1_binary, filtered_m2_binary)

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from matplotlib.legend_handler import HandlerLine2D
import matplotlib.patches as patches

#Enable LaTeX for all text in the plots
rc('text', usetex=True)
rc('font', family='serif')

# Define customization parameters
FONT_SIZE = 24        # Font size for labels and ticks
BORDER_THICKNESS = 1.5  # Thickness of the plot borders (spines)
LINE_WIDTH = 2.0       # Thickness of histogram lines
TICK_WIDTH = 1.5       # Thickness of the ticks

# Create subplots with a 2x2 layout
fig, ax = plt.subplots(2, 2, figsize=(12, 8))

# Plot the histograms in the subplots
bins = np.linspace(min(primary_mass), max(primary_mass), 20)
ax[0, 0].hist(primary_mass, bins=bins, histtype='step', color='green', alpha=0.5, linewidth=LINE_WIDTH)
ax[0, 0].hist(primary_mass_binary, bins=bins, histtype='step', color='purple', alpha=0.5, linewidth=LINE_WIDTH)
ax[0, 0].hist(filtered_primary_mass, bins=bins,  color='green', alpha=0.5, linewidth=LINE_WIDTH)
ax[0, 0].hist(filtered_primary_mass_binary, bins=bins,  color='purple', alpha=0.5, linewidth=LINE_WIDTH)
ax[0, 0].set_yscale('log')
ax[0, 0].set_xlabel(r'$\mathrm{M_1}$', fontsize=FONT_SIZE)
ax[0, 0].set_ylabel('N', fontsize=0.9*FONT_SIZE)

bins = np.linspace(min(mchirp), max(mchirp), 20)
ax[0, 1].hist(mchirp, bins=bins, histtype='step', color='green', alpha=0.5, linewidth=LINE_WIDTH)
ax[0, 1].hist(mchirp_binary, bins=bins, histtype='step', color='purple', alpha=0.5, linewidth=LINE_WIDTH)
ax[0, 1].hist(filtered_mchirp, bins=bins,  color='green', alpha=0.5, linewidth=LINE_WIDTH)
ax[0, 1].hist(filtered_mchirp_binary, bins=bins,  color='purple', alpha=0.5, linewidth=LINE_WIDTH)
ax[0, 1].set_yscale('log')
ax[0, 1].set_xlabel(r'$\mathrm{M_c}$', fontsize=0.9*FONT_SIZE)
ax[0, 1].set_ylabel('N', fontsize=FONT_SIZE)

bins = np.linspace(min(np.log10(gw_f)), max(np.log10(gw_f)), 20)
ax[1, 0].hist(np.log10(gw_f), bins=bins, histtype='step', color='green', alpha=0.5, linewidth=LINE_WIDTH)
ax[1, 0].hist(np.log10(gw_f_binary), bins=bins, histtype='step', color='purple', alpha=0.5, linewidth=LINE_WIDTH)
ax[1, 0].hist(np.log10(filtered_gw_f), bins=bins,  color='green', alpha=0.5, linewidth=LINE_WIDTH)
ax[1, 0].hist(np.log10(filtered_gw_f_binary), bins=bins,  color='purple', alpha=0.5, linewidth=LINE_WIDTH)
ax[1, 0].set_yscale('log')
ax[1, 0].set_xlabel(r'$\log_{10}(\mathrm{f_{gw}})$', fontsize=0.9*FONT_SIZE)
ax[1, 0].set_ylabel('N', fontsize=FONT_SIZE)

bins = np.linspace(min(np.log10(1 - e)), max(np.log10(1 - e)), 20)
ax[1, 1].hist(np.log10(1 - e), bins=bins, histtype='step', color='green', alpha=0.5, linewidth=LINE_WIDTH)
ax[1, 1].hist(np.log10(1 - e_binary), bins=bins, histtype='step', color='purple', alpha=0.5, linewidth=LINE_WIDTH)
ax[1, 1].hist(np.log10(1 - filtered_e), bins=bins,  color='green', alpha=0.5, linewidth=LINE_WIDTH)
ax[1, 1].hist(np.log10(1 - filtered_e_binary), bins=bins,  color='purple', alpha=0.5, linewidth=LINE_WIDTH)
ax[1, 1].set_yscale('log')
ax[1, 1].set_xlabel(r'$\log_{10}(\mathrm{1-e})$', fontsize=0.9*FONT_SIZE)
ax[1, 1].set_ylabel('N', fontsize=FONT_SIZE)

# Adjust border thickness and tick size
for a in ax.flat:
    # Customize borders (spines)
    for spine in a.spines.values():
        spine.set_linewidth(BORDER_THICKNESS)
    # Customize ticks
    a.tick_params(axis='both', which='major', labelsize=FONT_SIZE - 2, width=TICK_WIDTH)
    a.tick_params(axis='both', which='minor', labelsize=FONT_SIZE - 4, width=TICK_WIDTH)
    
legend_elements = [
    plt.Line2D([0], [0], color='green', lw=2 * LINE_WIDTH, label='triple'),
    plt.Line2D([0], [0], color='purple', lw=2 * LINE_WIDTH, label='binary'),
    patches.Patch(facecolor='green', edgecolor='black', label=r'$\mathrm{triple\,(\rho > 7)}$'),
    patches.Patch(facecolor='purple', edgecolor='black', label=r'$\mathrm{binary\,(\rho > 7)}$'),
]
# Set the position of the common legend
fig.legend(handles=legend_elements, fontsize=FONT_SIZE, handler_map={plt.Line2D: HandlerLine2D(numpoints=1)}, ncol=4, loc='upper center', bbox_to_anchor=(0.49, 1.0))# Adjust margins to prevent labels from being cut off
plt.subplots_adjust(top=0.9, bottom=0.1, hspace=0.3, wspace=0.3)
plt.savefig('pop_properties.pdf')

# Show the plot
plt.show()




