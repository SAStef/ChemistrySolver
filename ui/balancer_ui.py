"""
Terminal User Interface for Chemical Equation Balancer
"""
from ui.terminal_ui import display_title, display_results_header, wait_for_user
from chemistry_solver.balancer import parse_chemical_equation, balance_equation, format_balanced_equation, count_elements
from chemistry_solver.molar_mass import calculate_molar_mass

class BalancerUI:
    """UI class for chemical equation balancing."""
    
    def __init__(self):
        self.title = "CHEMICAL EQUATION BALANCER"
    
    def run(self):
        """Run the chemical equation balancer UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-3): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_simple_balance()
            elif choice == "2":
                self._handle_detailed_balance()
            elif choice == "3":
                self._handle_stoichiometry_calculation()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the equation balancer module menu."""
        menu = """
        [1] Balance chemical equation
        [2] Detailed equation balancing with element analysis
        [3] Stoichiometry calculations
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_simple_balance(self):
        """Handle basic equation balancing."""
        print("\n===== CHEMICAL EQUATION BALANCER =====")
        
        try:
            eq = input("Enter equation (e.g., H2 + O2 -> H2O): ")
            
            # Parse and balance the equation
            reactants, products = parse_chemical_equation(eq)
            balanced_reactants, balanced_products = balance_equation(reactants, products)
            balanced_equation = format_balanced_equation(balanced_reactants, balanced_products)
            
            display_results_header()
            print(f"Original equation: {eq}")
            print(f"Balanced equation: {balanced_equation}")
            
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_detailed_balance(self):
        """Handle equation balancing with detailed element analysis."""
        print("\n===== DETAILED EQUATION BALANCING =====")
        
        try:
            eq = input("Enter equation (e.g., H2 + O2 -> H2O): ")
            
            # Parse and balance the equation
            reactants, products = parse_chemical_equation(eq)
            balanced_reactants, balanced_products = balance_equation(reactants, products)
            balanced_equation = format_balanced_equation(balanced_reactants, balanced_products)
            
            display_results_header()
            print(f"Original equation: {eq}")
            print(f"Balanced equation: {balanced_equation}")
            
            # Element analysis - check conservation of mass
            print("\nElement Conservation Analysis:")
            print("-" * 50)
            print(f"{'Element':<10} {'Reactants':<15} {'Products':<15} {'Balanced?':<10}")
            print("-" * 50)
            
            # Count elements on both sides
            reactant_elements = {}
            for coef, formula in balanced_reactants:
                for element, count in count_elements(formula).items():
                    reactant_elements[element] = reactant_elements.get(element, 0) + coef * count
            
            product_elements = {}
            for coef, formula in balanced_products:
                for element, count in count_elements(formula).items():
                    product_elements[element] = product_elements.get(element, 0) + coef * count
            
            # Display element counts
            all_elements = sorted(set(list(reactant_elements.keys()) + list(product_elements.keys())))
            for element in all_elements:
                r_count = reactant_elements.get(element, 0)
                p_count = product_elements.get(element, 0)
                balanced = "✓" if r_count == p_count else "✗"
                print(f"{element:<10} {r_count:<15} {p_count:<15} {balanced:<10}")
            
            # Molecular weight analysis
            print("\nMolecular Weight Analysis:")
            print("-" * 50)
            
            # Calculate molecular weights
            reactant_weights = []
            total_reactant_weight = 0
            for coef, formula in balanced_reactants:
                result = calculate_molar_mass(formula)
                if result['success']:
                    weight = result['molar_mass']
                    reactant_weights.append((formula, coef, weight, coef * weight))
                    total_reactant_weight += coef * weight
            
            product_weights = []
            total_product_weight = 0
            for coef, formula in balanced_products:
                result = calculate_molar_mass(formula)
                if result['success']:
                    weight = result['molar_mass']
                    product_weights.append((formula, coef, weight, coef * weight))
                    total_product_weight += coef * weight
            
            # Display weights
            print("Reactants:")
            for formula, coef, weight, total in reactant_weights:
                print(f"  {formula}: {coef} × {weight:.4f} g/mol = {total:.4f} g/mol")
            print(f"Total reactant mass: {total_reactant_weight:.4f} g/mol")
            
            print("\nProducts:")
            for formula, coef, weight, total in product_weights:
                print(f"  {formula}: {coef} × {weight:.4f} g/mol = {total:.4f} g/mol")
            print(f"Total product mass: {total_product_weight:.4f} g/mol")
            
            # Check conservation of mass
            mass_difference = abs(total_reactant_weight - total_product_weight)
            if mass_difference < 0.001:
                print("\nMass is conserved! ✓")
            else:
                print(f"\nWarning: Mass difference of {mass_difference:.4f} g/mol detected")
            
        except Exception as e:
            print(f"Error: {str(e)}")
    
    def _handle_stoichiometry_calculation(self):
        """Handle stoichiometry calculations."""
        print("\n===== STOICHIOMETRY CALCULATOR =====")
        
        try:
            eq = input("Enter equation (e.g., H2 + O2 -> H2O): ")
            
            # Parse and balance the equation
            reactants, products = parse_chemical_equation(eq)
            balanced_reactants, balanced_products = balance_equation(reactants, products)
            balanced_equation = format_balanced_equation(balanced_reactants, balanced_products)
            
            display_results_header()
            print(f"Balanced equation: {balanced_equation}")
            
            # Get the known substance and amount
            print("\nEnter the known substance and amount:")
            
            all_substances = []
            for coef, formula in balanced_reactants:
                all_substances.append((formula, coef, "reactant"))
            for coef, formula in balanced_products:
                all_substances.append((formula, coef, "product"))
            
            print("\nAvailable substances:")
            for i, (formula, coef, side) in enumerate(all_substances, 1):
                print(f"  [{i}] {formula} (coefficient: {coef}, {side})")
            
            known_idx = int(input("\nSelect the number of the known substance: ")) - 1
            if known_idx < 0 or known_idx >= len(all_substances):
                raise ValueError("Invalid selection")
            
            known_formula, known_coef, _ = all_substances[known_idx]
            
            unit_type = input("Enter the unit type (moles/grams): ").lower()
            known_amount = float(input(f"Enter the amount of {known_formula} in {unit_type}: "))
            
            # Convert to moles if necessary
            known_moles = known_amount
            if unit_type == "grams":
                result = calculate_molar_mass(known_formula)
                if not result['success']:
                    raise ValueError(f"Could not calculate molar mass of {known_formula}")
                known_molar_mass = result['molar_mass']
                known_moles = known_amount / known_molar_mass
                print(f"\nConverting {known_amount} g of {known_formula} to moles:")
                print(f"Moles = grams / molar mass = {known_amount} / {known_molar_mass:.4f} = {known_moles:.6f} mol")
            
            # Get the target substance
            print("\nSelect target substance to calculate:")
            for i, (formula, coef, side) in enumerate(all_substances, 1):
                if i - 1 != known_idx:  # Skip the known substance
                    print(f"  [{i}] {formula} (coefficient: {coef}, {side})")
            
            target_idx = int(input("\nSelect the number of the target substance: ")) - 1
            if target_idx < 0 or target_idx >= len(all_substances) or target_idx == known_idx:
                raise ValueError("Invalid selection")
            
            target_formula, target_coef, _ = all_substances[target_idx]
            
            # Calculate mole ratio
            mole_ratio = target_coef / known_coef
            target_moles = known_moles * mole_ratio
            
            print(f"\nCalculating amount of {target_formula}:")
            print(f"Mole ratio: {target_coef} / {known_coef} = {mole_ratio}")
            print(f"Target moles = Known moles × Mole ratio = {known_moles:.6f} × {mole_ratio} = {target_moles:.6f} mol")
            
            # Calculate grams if requested
            output_unit = input("\nOutput unit (moles/grams): ").lower()
            if output_unit == "grams":
                result = calculate_molar_mass(target_formula)
                if not result['success']:
                    raise ValueError(f"Could not calculate molar mass of {target_formula}")
                target_molar_mass = result['molar_mass']
                target_grams = target_moles * target_molar_mass
                print(f"\nConverting {target_moles:.6f} mol of {target_formula} to grams:")
                print(f"Grams = moles × molar mass = {target_moles:.6f} × {target_molar_mass:.4f} = {target_grams:.4f} g")
                print(f"\nResult: {target_grams:.4f} g of {target_formula}")
            else:
                print(f"\nResult: {target_moles:.6f} mol of {target_formula}")
            
        except Exception as e:
            print(f"Error: {str(e)}")