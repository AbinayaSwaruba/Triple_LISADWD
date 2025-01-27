import numpy as np

# Define the parameters
Nb_in_range = 100000  # Example value for N_t, in range
fb_in_range = 0.1  # Example value for f_t, in range
fb = 0.5            # Example value for f_b
mt = 3.5           # Mass m_t (example value)
mb = 0.9            # Mass m_b (example value)
ms = 0.5            # Mass m_s (example value)

# Calculate Mtot, MSE, only binaries

Mtot_MSE = (Nb_in_range / (fb_in_range * fb)) * (fb * mb + (1 - fb) * ms)
# Display the result
print(f'Mtot, MSE: {Mtot_MSE}')

