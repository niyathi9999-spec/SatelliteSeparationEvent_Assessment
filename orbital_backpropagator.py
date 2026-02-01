"""The following program implements a numerical "back-propagator" 
using Cowell’s Method. It takes the current position and velocity
of several fragments and integrates their trajectories backward 
in time to find the exact moment and location where they all
"converge" to a single point. """

import numpy as np

"""This code uses scipy.integrate.solve_ivp to handle
the differential equations of motion in reverse."""

from scipy.integrate import solve_ivp

# Physical Constants
MU = 398600.4418  # Earth's gravitational parameter (km^3/s^2)

def equations_of_motion(t, state):
    """
    Defines the two-body differential equations: r'' = -mu * r / |r|^3
    'state' contains [x, y, z, vx, vy, vz]
    """
    r = state[:3]
    v = state[3:]
    r_mag = np.linalg.norm(r)
    
    # Acceleration vector
    a = -MU * r / (r_mag**3)
    
    # Return derivative [vx, vy, vz, ax, ay, az]
    return np.concatenate([v, a])

def backpropagate_fragments(fragment_states, rewind_time_seconds):
    """
    Rewinds fragment trajectories to find the convergence point.
    """
    history = {}
    
    # We integrate backward (from 0 to -rewind_time)
    t_span = (0, -rewind_time_seconds)
    t_eval = np.linspace(0, -rewind_time_seconds, 1000)
    
    for i, state in enumerate(fragment_states):
        sol = solve_ivp(equations_of_motion, t_span, state, t_eval=t_eval, rtol=1e-9)
        history[f"Frag_{i}"] = sol.y  # All positions/velocities over the rewind period
        
    return t_eval, history

# --- MOCK SCENARIO ---
# Current state of 3 fragments [x, y, z, vx, vy, vz] in km and km/s
# These were slightly separated by a small force 5 minutes (300s) ago
fragments = [
    [7000.1, 10.0, 5.0, 0.0, 7.501, 0.001],
    [6999.9, 10.1, 4.9, 0.001, 7.499, -0.002],
    [7000.0, 9.9, 5.1, -0.001, 7.500, 0.001]
]

rewind_duration = 600  # Look back 10 minutes
time_steps, results = backpropagate_fragments(fragments, rewind_duration)

# To find convergence: we look for the time 't' where the fragments are closest together.
# We'll calculate the mean distance from the centroid (dispersion).
min_dispersion = float('inf')
event_time_index = 0

for t_idx in range(len(time_steps)):
    # Get positions of all fragments at time t_idx
    # shape (3, N) -> 3 coordinates, N fragments
    pos_at_t = np.array([results[f][0:3, t_idx] for f in results])
    
    # Calculate centroid at this time
    centroid = np.mean(pos_at_t, axis=0)
    
    # Calculate distances of each fragment from the centroid
    distances = np.linalg.norm(pos_at_t - centroid, axis=1)
    
    # Measure of spread: mean distance from centroid
    dispersion = np.mean(distances)
    
    if dispersion < min_dispersion:
        min_dispersion = dispersion
        event_time_index = t_idx

# --- RESULTS ---
event_time = abs(time_steps[event_time_index])
# Centroid at the event time is the best estimate of location
event_loc = np.mean([results[f][0:3, event_time_index] for f in results], axis=0)

print(f"--- Separation Event Detection Result ---")
print(f"Estimated Event Time: {event_time:.1f} seconds ago")
print(f"Estimated ECI Location: {event_loc}")
print(f"Convergence Accuracy (Dispersion): {min_dispersion * 1000:.2f} meters")

"""How this Program Reaches the Solution
Numerical Integration (Cowell's Method): The solve_ivp function
acts as our time machine. By providing a negative time span, 
it calculates where the fragments must have been based on the laws of gravity.
State Vector Derivative: The equations_of_motion function ensures that
at every millisecond of the "rewind," the velocity and gravitational
acceleration are correctly updated. This handles Requirement 24’s need 
for "high-fidelity propagation."
Statistical Convergence Logic: The program doesn't just stop at a random time.
It iterates through the history and calculates the standard deviation of positions.
When the fragments are at their closest (the minimum dispersion point), 
that is mathematically defined as the Time of Separation.

Spatial Accuracy: The output gives the ECI (Earth-Centered Inertial) coordinates,
allowing the Product Architecture to cross-reference this location with the
Satellite Catalog to identify the "Parent" satellite."""