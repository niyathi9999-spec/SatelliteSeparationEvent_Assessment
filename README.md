# SatelliteSeparationEvent_Assessment
Code to submit for Digantara Assessment

The main back-propagation logic and example scenario: orbital_backpropagator.py
A simple separation-energy classifier example: separation-event.py

# orbital_backpropagator readme

## Key symbols (API)

orbital_backpropagator.equations_of_motion — computes derivatives for the two-body problem (state -> [v, a]).
orbital_backpropagator.backpropagate_fragments — integrates fragment states backward in time and returns time history.
orbital_backpropagator.MU — gravitational parameter used (km^3/s^2).
separation_event.calculate_separation_energy — computes fragment ΔV metrics and a basic classification.
How it works (concise)

The script integrates each fragment state backward in time using SciPy’s solve_ivp.
For each integration time-step it computes the centroid of fragment positions and the mean distance (dispersion).
The time of minimum dispersion is reported as the estimated separation time; the centroid at that time is the estimated ECI location.

## Outputs

Estimated event time (seconds ago)
Estimated ECI coordinates (km)
Convergence accuracy reported as mean dispersion (converted to meters in the printout)

# separation_event readme

## Purpose
Estimate how energetic a fragmentation/separation was by comparing fragment velocities to a chosen reference (centroid or parent-state). Useful as a quick classifier to distinguish low-energy detachments (e.g., deployment) from high-energy breakups.

## Inputs
- fragment_states: list or array of fragment state vectors.
  - Expected shape: (N, 6) or list of length-N where each element is [x, y, z, vx, vy, vz].
  - Units (recommended): positions in meters (m) or kilometers (km) and velocities in meters/second (m/s) or kilometers/second (km/s). Document and convert consistently before calling.
- reference_velocity (optional): 1x3 array used as the baseline velocity. If omitted the function uses the velocity centroid of the fragments.
- units (optional): string flag to indicate "m/s" or "km/s" so outputs are labeled correctly.
- thresholds (optional): dict with numeric thresholds for class labels (e.g., low/medium/high ΔV in chosen velocity units).

## Outputs
Returns a dictionary with:
- per_fragment_dv: array of ΔV magnitudes for each fragment (same units as input velocities).
- mean_dv: mean ΔV across fragments.
- max_dv: maximum ΔV across fragments.
- reference_velocity: velocity used as baseline (3-vector).
- classification: simple label (e.g., "low", "medium", "high") based on thresholds.
- metadata: units, number_of_fragments, timestamp (optional).

Printed output (if run as script) typically includes:
- Mean ΔV, Max ΔV
- Per-fragment ΔV list
- Classification summary

