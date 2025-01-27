import numpy as np

# Define the parameters
Nt_in_range = 100000  # Example value for N_t, in range
ft_in_range = 0.16   # Example value for f_t, in range
# Define the parameters
ft = 0.5            # Example value for f_t
fb = 0.3            # Example value for f_b
mt = 3.5           # Mass m_t (example value)
mb = 0.9            # Mass m_b (example value)
ms = 0.5            # Mass m_s (example value)

# Calculate Mtot, MSE
Mtot_MSE = (Nt_in_range / (ft_in_range * ft)) * (ft * mt + fb * mb + (1 - ft - fb) * ms)

# Display the result
print(f'Mtot, MSE: {Mtot_MSE}')
