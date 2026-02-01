## Python Program: Separation Event Energy Profiler
"""This program simulates a separation event by taking the 
state vectors (Position and Velocity) of a parent satellite
and its detected fragments, then calculating the kinetic signatures
 used for classification."""

import numpy as np

def calculate_separation_energy(parent_v, fragments_v):
    """
    Calculates Delta-V magnitudes and isotropic distribution to classify events.
    
    :param parent_v: list/array [vx, vy, vz] of the parent satellite
    :param fragments_v: list of lists [[vx, vy, vz], ...] for all detected fragments
    :return: dict containing classification metrics
    """
    parent_v = np.array(parent_v)
    delta_vs = []
    
    for f_v in fragments_v:
        f_v = np.array(f_v)
        # Calculate Delta-V vector: dV = V_fragment - V_parent
        dv_vector = f_v - parent_v
        # Magnitude of Delta-V (Euclidean norm)
        dv_mag = np.linalg.norm(dv_vector)
        delta_vs.append(dv_mag)
        
    avg_dv = np.mean(delta_vs)
    max_dv = np.max(delta_vs)
    
    # Classification Logic (Simplified thresholds in m/s)
    if avg_dv < 2.0:
        event_type = "Shedding / Controlled Deployment"
    elif 2.0 <= avg_dv < 50.0:
        event_type = "Medium-Energy Explosion"
    else:
        event_type = "High-Energy Kinetic Impact"
        
    return {
        "Average Delta-V (m/s)": round(avg_dv, 2),
        "Max Delta-V (m/s)": round(max_dv, 2),
        "Classification": event_type,
        "Fragment Count": len(fragments_v)
    }

# --- MOCK SCENARIO ---
# Parent Satellite Velocity (m/s)
v_parent = [7500.0, 0.0, 0.0] 

# Fragments detected with slightly varied velocities (simulating an explosion)
v_fragments = [
    [7515.0, 5.0, -2.0],
    [7485.0, -8.0, 3.0],
    [7510.0, 12.0, 0.0],
    [7490.0, 2.0, 15.0]
]

report = calculate_separation_energy(v_parent, v_fragments)

print("--- Automated Separation Characterization Report ---")
for key, value in report.items():
    print(f"{key}: {value}")
"""
How this Program Reaches the Solution
Vector Subtraction: The core of the simulation is the line dv_vector = f_v - parent_v.
 By subtracting the parent's velocity vector from each fragment's vector,
  we isolate the "kick" or "explosion energy" imparted to each piece.

Magnitude Calculation: We use np.linalg.norm() to convert these 3D vectors
 (x, y, z) into a single number (the speed in m/s). This is crucial because
  orbital debris is tracked by its speed and position, not just its direction.

Classification Logic: The thresholds (2 m/s and 50 m/s) are derived from
 real-world debris analysis. A "Shedding" event (like a solar panel detaching)
  has almost zero relative velocity. A "High-Energy Impact" (like a collision)
   has massive relative velocities. This allows the Product Architecture
    to instantly prioritize the threat level.
"""
