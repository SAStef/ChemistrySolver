from molar_mass import calculate_molar_mass
from balancer import parse_chemical_equation, balance_equation, format_balanced_equation
from stoichiometry import (solve_stoichiometry_problem, solve_multireactant_problem, 
                        solve_gas_stoichiometry_problem)
from acid_base import identify_acid_base, analyze_compound_list
from oxidation_state import calculate_oxidation_number, display_oxidation_result
from redox_reactions import (balance_redox_reaction, identify_oxidation_changes, 
                            find_molar_ratio,
                            determine_reaction_favorability)

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
        print("0. Exit")
        choice = input("Enter choice: ")
        action = actions.get(choice)
        if action:
            action()
        else:
            print("Invalid choice.")

if __name__ == "__main__":
    main()