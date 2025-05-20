from molar_mass import calculate_molar_mass
from balancer import parse_chemical_equation, balance_equation, format_balanced_equation
from stoichiometry import solve_stoichiometry_problem, solve_multireactant_problem

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

def main():
    actions = {
        '1': handle_molar_mass,
        '2': handle_balance,
        '3': handle_stoichiometry,
        '4': handle_multireactant,
        '0': exit
    }
    while True:
        print("\n===== Chemical Equation Solver =====")
        print("1. Calculate molar mass")
        print("2. Balance a chemical equation")
        print("3. Solve stoichiometry problem")
        print("4. Solve multireactant problem (limiting reactant)")
        print("0. Exit")
        choice = input("Enter choice: ")
        action = actions.get(choice)
        if action:
            action()
        else:
            print("Invalid choice.")

if __name__ == "__main__":
    main()