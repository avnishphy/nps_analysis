import math

def calculate_kinematics(beam_energy, proton_angle):
    """
    Calculate the scattered proton momentum, scattered electron momentum, and scattered electron angle.

    Parameters:
    beam_energy (float): Initial beam energy in GeV.
    proton_angle (float): Scattered proton angle (theta_p) in degrees.

    Returns:
    tuple: Scattered proton momentum (P_p), scattered electron momentum (P_e), and scattered electron angle (theta_e) in degrees.
    """
    mass_e = 0.000511  # Electron mass in GeV
    mass_p = 0.938272  # Proton mass in GeV

    # Initial electron momentum
    initial_momentum = math.sqrt(beam_energy**2 - mass_e**2)

    # Convert proton angle from degrees to radians
    theta_p_rad = math.radians(proton_angle)

    # Compute the final proton energy correctly
    E_p_final = beam_energy + mass_p - math.sqrt(initial_momentum**2 + mass_e**2)

    # Ensure E_p_final > mass_p to avoid sqrt of negative number
    if E_p_final < mass_p:
        raise ValueError("Invalid energy-momentum balance: Check input values.")

    # Compute scattered proton momentum
    P_p = math.sqrt(E_p_final**2 - mass_p**2)

    # Compute scattered electron momentum using energy conservation
    E_e_final = beam_energy + mass_p - E_p_final
    P_e = math.sqrt(E_e_final**2 - mass_e**2)

    # Compute scattered electron angle using momentum conservation
    sin_theta_e = (P_p * math.sin(theta_p_rad)) / P_e

    # Ensure valid range for asin
    if abs(sin_theta_e) > 1:
        raise ValueError("Invalid input: sin(theta_e) out of range. Check kinematics.")

    theta_e_rad = math.asin(sin_theta_e)
    theta_e_deg = math.degrees(theta_e_rad)

    return P_p, P_e, theta_e_deg

def main():
    # Get dynamic input from the user
    beam_energy = float(input("Enter the initial beam energy in GeV: "))
    proton_angle = float(input("Enter the scattered proton angle (theta_p) in degrees: "))

    try:
        # Calculate kinematic variables
        proton_momentum, electron_momentum, electron_angle = calculate_kinematics(beam_energy, proton_angle)

        # Output the results
        print(f"\nScattered Proton Momentum (P_p): {proton_momentum:.6f} GeV/c")
        print(f"Scattered Electron Momentum (P_e): {electron_momentum:.6f} GeV/c")
        print(f"Scattered Electron Angle (theta_e): {electron_angle:.2f} degrees")

    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
