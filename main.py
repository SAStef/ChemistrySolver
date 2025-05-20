from molar_mass import calculate_molar_mass
from balancer import parse_chemical_equation, balance_equation, format_balanced_equation
from stoichiometry import (solve_stoichiometry_problem, solve_multireactant_problem, 
                        solve_gas_stoichiometry_problem)
from acid_base import identify_acid_base, analyze_compound_list
from oxidation_state import calculate_oxidation_number, display_oxidation_result
from redox_reactions import (balance_redox_reaction, identify_oxidation_changes, 
                            find_molar_ratio,
                            determine_reaction_favorability)
from functional_groups import (identify_functional_groups, check_functional_groups,
                            find_missing_functional_groups, explain_functional_groups_in_compound,
                            solve_functional_group_problem)
# Import thermodynamics module
from thermodynamics import (calculate_heat, calculate_temperature_change, calculate_final_temperature,
                           calculate_molar_heat, solve_thermal_equilibrium, solve_thermal_equilibrium_with_molar_heat,
                           handle_heat_transfer_problem, handle_heat_transfer_with_molar_heat, solve_mixture_problem)

from colligative_properties import (
    solve_freezing_point_problem, calculate_boiling_point_elevation,
    calculate_osmotic_pressure, calculate_vapor_pressure_lowering,
    print_solution, FREEZING_POINT_CONSTANTS, BOILING_POINT_CONSTANTS
)

def handle_freezing_point_depression():
    """
    Handler function for solving molecular weight problems using freezing point depression.
    """
    print("\n=== Freezing Point Depression Calculator ===")
    solvent = input("Enter solvent formula (e.g., 'H2O', 'C6H6' for benzene): ")
    
    # Try to lookup Kf in constants
    k_f = None
    solvent_lower = solvent.lower()
    for key, value in FREEZING_POINT_CONSTANTS.items():
        if key in solvent_lower or solvent_lower in key:
            k_f = value
            print(f"Found freezing point constant for {key}: {value} °C/m")
            use_constant = input("Use this constant? (y/n): ").lower() == 'y'
            if use_constant:
                break
            else:
                k_f = None
    
    if k_f is None:
        k_f = float(input("Enter freezing point depression constant (Kf in °C/m): "))
    
    T_pure = float(input("Enter freezing point of pure solvent (°C): "))
    T_solution = float(input("Enter freezing point of solution (°C): "))
    solute_mass = float(input("Enter mass of solute (g): "))
    solvent_mass = float(input("Enter mass of solvent (g): "))
    
    ion_input = input("Does the solute ionize? (y/n): ").lower()
    ionization_factor = 1
    if ion_input == 'y':
        ionization_factor = float(input("Enter van 't Hoff factor (i): "))
    
    result = solve_freezing_point_problem(
        T_pure=T_pure,
        T_solution=T_solution,
        solvent=solvent,
        K_f=k_f,
        solute_mass=solute_mass,
        solvent_mass=solvent_mass,
        ionization_factor=ionization_factor
    )
    
    print_solution(result, "Freezing Point Depression")

def handle_boiling_point_elevation():
    """
    Handler function for solving molecular weight problems using boiling point elevation.
    """
    print("\n=== Boiling Point Elevation Calculator ===")
    solvent = input("Enter solvent formula (e.g., 'H2O', 'C6H6' for benzene): ")
    
    # Try to lookup Kb in constants
    k_b = None
    solvent_lower = solvent.lower()
    for key, value in BOILING_POINT_CONSTANTS.items():
        if key in solvent_lower or solvent_lower in key:
            k_b = value
            print(f"Found boiling point constant for {key}: {value} °C/m")
            use_constant = input("Use this constant? (y/n): ").lower() == 'y'
            if use_constant:
                break
            else:
                k_b = None
    
    if k_b is None:
        k_b = float(input("Enter boiling point elevation constant (Kb in °C/m): "))
    
    T_pure = float(input("Enter boiling point of pure solvent (°C): "))
    T_solution = float(input("Enter boiling point of solution (°C): "))
    solute_mass = float(input("Enter mass of solute (g): "))
    solvent_mass = float(input("Enter mass of solvent (g): "))
    
    ion_input = input("Does the solute ionize? (y/n): ").lower()
    ionization_factor = 1
    if ion_input == 'y':
        ionization_factor = float(input("Enter van 't Hoff factor (i): "))
    
    result = calculate_boiling_point_elevation(
        T_pure=T_pure,
        T_solution=T_solution,
        K_b=k_b,
        solute_mass=solute_mass,
        solvent=solvent,
        solvent_mass=solvent_mass,
        ionization_factor=ionization_factor
    )
    
    print_solution(result, "Boiling Point Elevation")

def handle_osmotic_pressure():
    """
    Handler function for solving molecular weight problems using osmotic pressure.
    """
    print("\n=== Osmotic Pressure Calculator ===")
    
    osmotic_pressure_atm = float(input("Enter osmotic pressure (atm): "))
    temperature_c = float(input("Enter temperature (°C): "))
    solution_volume_L = float(input("Enter solution volume (L): "))
    solute_mass = float(input("Enter mass of solute (g): "))
    
    ion_input = input("Does the solute ionize? (y/n): ").lower()
    ionization_factor = 1
    if ion_input == 'y':
        ionization_factor = float(input("Enter van 't Hoff factor (i): "))
    
    result = calculate_osmotic_pressure(
        osmotic_pressure_atm=osmotic_pressure_atm,
        temperature_c=temperature_c,
        solution_volume_L=solution_volume_L,
        solute_mass=solute_mass,
        ionization_factor=ionization_factor
    )
    
    print_solution(result, "Osmotic Pressure")

def handle_vapor_pressure_lowering():
    """
    Handler function for solving molecular weight problems using vapor pressure lowering.
    """
    print("\n=== Vapor Pressure Lowering Calculator ===")
    
    solvent = input("Enter solvent formula (e.g., 'H2O', 'C6H6' for benzene): ")
    P_pure = float(input("Enter vapor pressure of pure solvent: "))
    P_solution = float(input("Enter vapor pressure of solution: "))
    solute_mass = float(input("Enter mass of solute (g): "))
    solvent_mass = float(input("Enter mass of solvent (g): "))
    
    result = calculate_vapor_pressure_lowering(
        P_pure=P_pure,
        P_solution=P_solution,
        solute_mass=solute_mass,
        solvent=solvent,
        solvent_mass=solvent_mass
    )
    
    print_solution(result, "Vapor Pressure Lowering")
    
def handle_molar_mass():
    formula = input("Enter chemical formula: ")
    result = calculate_molar_mass(formula)
    if result['success']:
        print(f"Molar mass: {result['molar_mass']:.4f} g/mol")
        for e in result['composition']:
            print(f"  {e['element']}: {e['count']} × {e['atomic_mass']:.4f} u = {e['contribution']:.4f}")
    else:
        print(f"Error: {result['error']}")

def handle_balance():
    eq = input("Enter equation (e.g., H2 + O2 -> H2O): ")
    try:
        r, p = parse_chemical_equation(eq)
        br, bp = balance_equation(r, p)
        print("Balanced equation:", format_balanced_equation(br, bp))
    except Exception as e:
        print("Error:", str(e))

def handle_stoichiometry():
    eq = input("Equation: ")
    given = input("Known compound: ")
    mass = float(input("Mass in g: "))
    target = input("Target compound: ")
    try:
        result = solve_stoichiometry_problem(eq, given, mass, target)
        print("\n".join(result["steps"]))
        print(f"\nMass of {target}: {result['target_mass']:.4f} g")
    except Exception as e:
        print("Error:", str(e))

def handle_gas_stoichiometry():
    eq = input("Equation: ")
    given = input("Known compound: ")
    mass = float(input("Mass in g: "))
    target = input("Target gas compound: ")
    temp = float(input("Temperature (°C): "))
    pressure = float(input("Pressure (atm): "))
    
    try:
        result = solve_gas_stoichiometry_problem(eq, given, mass, target, temp, pressure)
        print("\n".join(result["steps"]))
        print(f"\nVolume of {target} gas: {result['gas_volume']:.4f} L")
    except Exception as e:
        print("Error:", str(e))

def handle_multireactant():
    eq = input("Equation: ")
    reactant_count = int(input("Number of reactants: "))
    reactant_data = {}
    
    for i in range(reactant_count):
        reactant = input(f"Reactant {i+1} formula: ")
        mass = float(input(f"Mass of {reactant} in g: "))
        reactant_data[reactant] = mass
    
    target = input("Target compound: ")
    
    try:
        result = solve_multireactant_problem(eq, reactant_data, target)
        print("\n".join(result["steps"]))
        print(f"\nLimiting reactant: {result['limiting_reactant']}")
        print(f"Mass of {target}: {result['target_mass']:.4f} g")
    except Exception as e:
        print("Error:", str(e))

def handle_acid_base_single():
    compound = input("Enter a chemical formula (e.g., HCl, NaOH): ")
    result = identify_acid_base(compound)
    print(f"\nCompound: {compound}")
    print(f"Classification: {result['classification']}")
    print(f"Explanation: {result['explanation']}")

def handle_acid_base_list():
    compound_input = input("Enter compounds separated by commas (e.g., HCl, H2SO4, CH3COOH): ")
    compounds = [comp.strip() for comp in compound_input.split(",")]
    analyze_compound_list(compounds)
    
def handle_oxidation_number():
    compound = input("Enter a chemical compound (e.g., CrO2Cl2, Fe2O3): ")
    element = input("Enter the element to find oxidation number for: ")
    try:
        result = calculate_oxidation_number(compound, element)
        display_oxidation_result(result)
    except Exception as e:
        print(f"Error: {str(e)}")

def handle_redox_balance():
    equation = input("Enter the redox reaction equation: ")
    environment = input("Environment (acidic/basic): ").lower()
    
    if environment not in ["acidic", "basic"]:
        environment = "acidic"  # default to acidic
    
    try:
        result = balance_redox_reaction(equation, environment)
        
        print("\n=== Balanced Redox Reaction ===")
        print(f"Balanced equation: {result['balanced_equation']}")
        print(f"Environment: {result['environment']}")
        
        print("\nElements undergoing redox:")
        for element in result['redox_elements']:
            change = element['change']
            direction = "oxidized" if change > 0 else "reduced"
            print(f"  {element['element']}: {element['reactant_oxidation']} → {element['product_oxidation']} ({direction})")
        
        print("\nOxidizing agents:", ", ".join(result['oxidizing_agents']))
        print("Reducing agents:", ", ".join(result['reducing_agents']))
    except Exception as e:
        print(f"Error: {str(e)}")

def handle_reaction_favorability():
    equation = input("Enter the redox reaction equation: ")
    
    try:
        result = determine_reaction_favorability(equation)
        
        print("\n=== Reaction Favorability Analysis ===")
        print(f"Balanced equation: {result['balanced_equation']}")
        
        if result['cell_potential'] is not None:
            print(f"\nStandard Cell Potential (E°cell): {result['cell_potential']:.2f} V")
            
            if result['favorable']:
                print("\nResult: FAVORABLE at standard conditions")
            else:
                print("\nResult: NOT FAVORABLE at standard conditions")
            
            print("\nExplanation:")
            print(result['explanation'])
        else:
            print("\nCould not determine favorability. Missing standard potentials for one or more half-reactions.")
    except Exception as e:
        print(f"Error: {str(e)}")

def handle_molar_ratio():
    equation = input("Enter the balanced equation: ")
    compound1 = input("First compound: ")
    compound2 = input("Second compound: ")
    
    try:
        ratio = find_molar_ratio(equation, compound1, compound2)
        
        if ratio is not None:
            print(f"Molar ratio of {compound1} to {compound2}: {ratio:.2f}:{1}")
            print(f"For every {ratio:.2f} moles of {compound1}, you need 1 mole of {compound2}.")
        else:
            print("Could not determine molar ratio. Make sure both compounds are in the equation.")
    except Exception as e:
        print(f"Error: {str(e)}")

# New handlers for functional groups
def handle_identify_functional_groups():
    """
    Handler function for identifying functional groups in a compound.
    """
    print("\n=== Identify Functional Groups ===")
    compound = input("Enter compound name: ")
    result = explain_functional_groups_in_compound(compound_name=compound)
    
    print(f"\nCompound: {result['compound']}")
    if result['functional_groups']:
        print("Identified functional groups:")
        for i, (group, explanation) in enumerate(zip(result['functional_groups'], result['explanations'])):
            print(f"  {i+1}. {group.capitalize()}: {explanation}")
    else:
        print("No functional groups identified.")

def handle_functional_group_problem():
    """
    Handler function for solving which functional group is NOT present in a compound.
    """
    print("\n=== Functional Group Analysis ===")
    compound = input("Enter compound name: ")
    print("Enter functional groups to check (comma-separated):")
    options_input = input("e.g., methyl, carboxyl, hydroxyl, amine, halogen: ")
    options = [opt.strip().lower() for opt in options_input.split(",")]
    
    result = solve_functional_group_problem(compound, options)
    
    print(f"\n=== Functional Group Analysis for {result['compound']} ===")
    print("Present functional groups:")
    for group in result['present_groups']:
        print(f"  - {group}")
    
    print("\nMissing functional groups:")
    for group in result['missing_groups']:
        print(f"  - {group}")
    
    if len(result['missing_groups']) == 1:
        print(f"\nAnswer: {result['missing_groups'][0]} is NOT present in {result['compound']}.")
    elif len(result['missing_groups']) > 1:
        print(f"\nAnswer: Multiple functional groups are not present in {result['compound']}.")
    else:
        print(f"\nAnswer: All specified functional groups are present in {result['compound']}.")

def handle_check_specific_functional_groups():
    """
    Handler function for checking specific functional groups in a compound.
    """
    print("\n=== Check Specific Functional Groups ===")
    compound = input("Enter compound name: ")
    print("Enter specific functional groups to check (comma-separated):")
    groups_input = input("e.g., methyl, carboxyl, hydroxyl, amine, halogen: ")
    groups_to_check = [group.strip().lower() for group in groups_input.split(",")]
    
    result = check_functional_groups(compound_name=compound, groups_to_check=groups_to_check)
    
    print(f"\nResults for {compound}:")
    for group, present in result.items():
        status = "Present" if present else "Not present"
        print(f"  {group.capitalize()}: {status}")

# New handlers for thermodynamics
def handle_heat_calculation():
    """
    Handler function for calculating heat using q = m × c × ΔT.
    """
    print("\n=== Heat Transfer Calculation ===")
    mass = float(input("Enter mass (g): "))
    specific_heat = float(input("Enter specific heat capacity (J/(g·K)): "))
    initial_temp = float(input("Enter initial temperature (°C): "))
    final_temp = float(input("Enter final temperature (°C): "))
    
    delta_t = final_temp - initial_temp
    heat = calculate_heat(mass, specific_heat, delta_t)
    
    print(f"\nResults:")
    print(f"Temperature change (ΔT): {delta_t:.2f} °C")
    print(f"Heat transferred (q): {heat:.2f} J")
    print(f"Heat transferred: {heat/1000:.2f} kJ")

def handle_temperature_change_calculation():
    """
    Handler function for calculating temperature change using ΔT = q / (m × c).
    """
    print("\n=== Temperature Change Calculation ===")
    heat = float(input("Enter heat energy (J): "))
    mass = float(input("Enter mass (g): "))
    specific_heat = float(input("Enter specific heat capacity (J/(g·K)): "))
    
    delta_t = calculate_temperature_change(heat, mass, specific_heat)
    
    print(f"\nResults:")
    print(f"Temperature change (ΔT): {delta_t:.2f} °C")

def handle_thermal_equilibrium():
    """
    Handler function for solving thermal equilibrium problems.
    """
    print("\n=== Thermal Equilibrium Problem ===")
    print("This will calculate the final temperature when two substances reach thermal equilibrium.")
    
    # Substance 1
    print("\nSubstance 1:")
    mass1 = float(input("Enter mass (g): "))
    specific_heat1 = float(input("Enter specific heat capacity (J/(g·K)): "))
    initial_temp1 = float(input("Enter initial temperature (°C): "))
    
    # Substance 2
    print("\nSubstance 2:")
    mass2 = float(input("Enter mass (g): "))
    specific_heat2 = float(input("Enter specific heat capacity (J/(g·K)): "))
    initial_temp2 = float(input("Enter initial temperature (°C): "))
    
    # Solve problem
    result = handle_heat_transfer_problem(
        mass1, specific_heat1, initial_temp1,
        mass2, specific_heat2, initial_temp2
    )
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
    print(f"Heat transferred by substance 1: {result['heat_transferred_1']:.2f} J")
    print(f"Heat transferred by substance 2: {result['heat_transferred_2']:.2f} J")

def handle_thermal_equilibrium_molar():
    """
    Handler function for solving thermal equilibrium problems using molar heat capacities.
    """
    print("\n=== Thermal Equilibrium Problem (Molar Heat Capacities) ===")
    print("This will calculate the final temperature when two substances reach thermal equilibrium.")
    
    # Substance 1
    print("\nSubstance 1:")
    name1 = input("Enter substance name or chemical formula: ")
    mass1 = float(input("Enter mass (g): "))
    
    # Try to calculate molar mass if a chemical formula is given
    try:
        molar_mass_result = calculate_molar_mass(name1)
        if molar_mass_result['success']:
            molar_mass1 = molar_mass_result['molar_mass']
            print(f"Calculated molar mass: {molar_mass1:.4f} g/mol")
        else:
            print("Could not automatically calculate molar mass.")
            molar_mass1 = float(input("Enter molar mass (g/mol): "))
    except:
        print("Could not automatically calculate molar mass.")
        molar_mass1 = float(input("Enter molar mass (g/mol): "))
    
    molar_heat_capacity1 = float(input("Enter molar heat capacity (J/(mol·K)): "))
    initial_temp1 = float(input("Enter initial temperature (°C): "))
    
    # Substance 2
    print("\nSubstance 2:")
    name2 = input("Enter substance name or chemical formula: ")
    mass2 = float(input("Enter mass (g): "))
    
    # Try to calculate molar mass if a chemical formula is given
    try:
        molar_mass_result = calculate_molar_mass(name2)
        if molar_mass_result['success']:
            molar_mass2 = molar_mass_result['molar_mass']
            print(f"Calculated molar mass: {molar_mass2:.4f} g/mol")
        else:
            print("Could not automatically calculate molar mass.")
            molar_mass2 = float(input("Enter molar mass (g/mol): "))
    except:
        print("Could not automatically calculate molar mass.")
        molar_mass2 = float(input("Enter molar mass (g/mol): "))
    
    molar_heat_capacity2 = float(input("Enter molar heat capacity (J/(mol·K)): "))
    initial_temp2 = float(input("Enter initial temperature (°C): "))
    
    # Solve problem
    result = handle_heat_transfer_with_molar_heat(
        mass1, molar_mass1, molar_heat_capacity1, initial_temp1,
        mass2, molar_mass2, molar_heat_capacity2, initial_temp2
    )
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
    print(f"Heat transferred by substance 1: {result['heat_transferred_1']:.2f} J")
    print(f"Heat transferred by substance 2: {result['heat_transferred_2']:.2f} J")

def handle_mixture_thermal_equilibrium():
    """
    Handler function for solving thermal equilibrium problems with multiple substances.
    """
    print("\n=== Multiple Substance Thermal Equilibrium Problem ===")
    print("This will calculate the final temperature when multiple substances reach thermal equilibrium.")
    
    substance_count = int(input("Enter number of substances: "))
    substances = []
    
    for i in range(substance_count):
        print(f"\nSubstance {i+1}:")
        name = input("Enter substance name or chemical formula: ")
        mass = float(input("Enter mass (g): "))
        
        # Try to automatically calculate molar mass
        use_molar_heat = input("Use molar heat capacity? (y/n): ").lower() == 'y'
        
        if use_molar_heat:
            # Try to calculate molar mass if a chemical formula is given
            try:
                molar_mass_result = calculate_molar_mass(name)
                if molar_mass_result['success']:
                    molar_mass = molar_mass_result['molar_mass']
                    print(f"Calculated molar mass: {molar_mass:.4f} g/mol")
                else:
                    print("Could not automatically calculate molar mass.")
                    molar_mass = float(input("Enter molar mass (g/mol): "))
            except:
                print("Could not automatically calculate molar mass.")
                molar_mass = float(input("Enter molar mass (g/mol): "))
                
            molar_heat_capacity = float(input("Enter molar heat capacity (J/(mol·K)): "))
            # Calculate specific heat
            specific_heat = molar_heat_capacity / molar_mass
            print(f"Calculated specific heat: {specific_heat:.4f} J/(g·K)")
        else:
            specific_heat = float(input("Enter specific heat capacity (J/(g·K)): "))
        
        initial_temp = float(input("Enter initial temperature (°C): "))
        
        substance = {
            'name': name,
            'mass': mass,
            'specific_heat': specific_heat,
            'initial_temp': initial_temp
        }
        
        if use_molar_heat:
            substance['molar_mass'] = molar_mass
            substance['molar_heat_capacity'] = molar_heat_capacity
        
        substances.append(substance)
    
    # Solve problem
    result = solve_mixture_problem(substances)
    
    print("\n=== Solution ===")
    for step in result["steps"]:
        print(step)
    
    print(f"\nFinal equilibrium temperature: {result['final_temp']:.2f} °C")
    print("\nHeat transferred by each substance:")
    for transfer in result["heat_transfers"]:
        print(f"  {transfer['name']}: {transfer['heat']:.2f} J (ΔT = {transfer['delta_t']:.2f} °C)")

def main():
    actions = {
        '1': handle_molar_mass,
        '2': handle_balance,
        '3': handle_stoichiometry,
        '4': handle_multireactant,
        '5': handle_acid_base_single,
        '6': handle_acid_base_list,
        '7': handle_oxidation_number,
        '8': handle_redox_balance,
        '9': handle_molar_ratio,
        '10': handle_gas_stoichiometry,
        '11': handle_reaction_favorability,
        '12': handle_identify_functional_groups,
        '13': handle_check_specific_functional_groups,
        '14': handle_functional_group_problem,
        '15': handle_freezing_point_depression,  
        '16': handle_boiling_point_elevation,    
        '17': handle_osmotic_pressure,           
        '18': handle_vapor_pressure_lowering,    
        '0': exit
    }
    while True:
        print("\n===== Chemical Equation Solver =====")
        print("1. Calculate molar mass")
        print("2. Balance a chemical equation")
        print("3. Solve stoichiometry problem")
        print("4. Solve multireactant problem (limiting reactant)")
        print("5. Identify acid/base compound")
        print("6. Analyze list of compounds")
        print("7. Calculate oxidation number")
        print("8. Balance a redox reaction")
        print("9. Find molar ratio between compounds")
        print("10. Solve gas stoichiometry problem")
        print("11. Determine redox reaction favorability")
        print("12. Identify functional groups in a compound")
        print("13. Check specific functional groups")
        print("14. Find missing functional groups in a compound")
        print("---- Colligative Properties ----")
        print("15. Freezing point depression (molecular weight)")
        print("16. Boiling point elevation (molecular weight)")
        print("17. Osmotic pressure (molecular weight)")
        print("18. Vapor pressure lowering (molecular weight)")
        print("0. Exit")
        
        choice = input("Enter choice: ")
        action = actions.get(choice)
        if action:
            action()
        else:
            print("Invalid choice.")

if __name__ == "__main__":
    main()