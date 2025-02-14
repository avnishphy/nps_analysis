import math

def calculate_proton_properties(beam_energy, electron_angle, electron_momentum):
    """
    Calculate the proton angle and momentum based on the provided inputs.

    Parameters:
    beam_energy (float): Initial beam energy.
    electron_angle (float): Scattered electron angle (theta_e) in degrees.
    electron_momentum (float): Scattered electron momentum (P_e).

    Returns:
    tuple: Proton angle (theta_p) in degrees and proton momentum (P_p).
    """
    mass_e = 0.000511 # in GeV
    mass_p = 0.938272 # in GeV

    initial_momentum = math.sqrt(beam_energy**2 - mass_e**2)
    print("initial momentum: ", initial_momentum)

    # Convert the electron angle from degrees to radians
    theta_e = math.radians(electron_angle)

    calculated_scattered_elec_mom = initial_momentum/(1 + (2*initial_momentum/mass_p)*(math.sin(theta_e/2))**2)

    # Calculate the proton angle theta_p using the given formula
    tan_theta_p = (calculated_scattered_elec_mom * math.sin(theta_e)) / (initial_momentum - calculated_scattered_elec_mom * math.cos(theta_e))
    theta_p = math.atan(tan_theta_p)
    # theta_p = math.atan(1/((1 + (initial_momentum/mass_p))*math.tan(theta_e/2)))


    # Convert the proton angle from radians to degrees
    theta_p_degrees = math.degrees(theta_p)

    # Calculate the proton momentum P_p using the given formula
    # proton_momentum = (electron_momentum * math.sin(theta_e)) / math.sin(theta_p)

    # proton_momentum = math.sqrt((initial_momentum + mass_p - electron_momentum)**2 - mass_p**2)
    proton_momentum = math.sqrt((initial_momentum + mass_p - calculated_scattered_elec_mom)**2 - mass_p**2)

    return theta_p_degrees, proton_momentum, calculated_scattered_elec_mom

def main():
    # Get dynamic input from the user
    beam_energy = float(input("Enter the initial beam energy in GeV: "))
    electron_angle = float(input("Enter the scattered electron angle (theta_e) in degrees: "))
    electron_momentum = float(input("Enter the scattered electron momentum (P_e) in GeV/c: "))

    # Calculate the proton angle and momentum
    proton_angle, proton_momentum, calculated_scattered_elec_mom = calculate_proton_properties(beam_energy, electron_angle, electron_momentum)

    # Output the results
    print(f"\nProton Angle (theta_p): {proton_angle:.2f} degrees")
    print(f"Proton Momentum (P_p): {proton_momentum:.6f} GeV/c")
    print(f"Calculated scattered electron momentum: {calculated_scattered_elec_mom:.6f} GeV/c")

if __name__ == "__main__":
    main()
