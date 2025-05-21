"""
Terminal User Interface for Chemical Kinetics Problems
"""
from ui.terminal_ui import display_title, display_results_header, display_steps, wait_for_user
from chemistry_solver.kinetics import (
    solve_first_order_kinetics, 
    calculate_rate_constant_from_half_life,
    calculate_half_life_from_rate_constant, 
    calculate_concentration_after_time,
    calculate_fraction_remaining, 
    as_simplified_fraction
)

class KineticsUI:
    """UI class for chemical kinetics problems."""
    
    def __init__(self):
        self.title = "CHEMICAL KINETICS PROBLEM SOLVER"
    
    def run(self):
        """Run the chemical kinetics UI."""
        display_title(self.title)
        
        while True:
            self._display_menu()
            choice = input("\nEnter choice (0-5): ").strip()
            
            if choice == "0":
                # Return to main menu
                return
            elif choice == "1":
                self._handle_rate_constant_calculation()
            elif choice == "2":
                self._handle_half_life_calculation()
            elif choice == "3":
                self._handle_concentration_calculation()
            elif choice == "4":
                self._handle_time_calculation()
            elif choice == "5":
                self._handle_complete_analysis()
            else:
                print("Invalid choice. Please try again.")
            
            wait_for_user()
    
    def _display_menu(self):
        """Display the chemical kinetics module menu."""
        menu = """
        [1] Calculate rate constant from half-life
        [2] Calculate half-life from rate constant
        [3] Calculate concentration after time
        [4] Calculate time to reach concentration
        [5] Complete first-order kinetics analysis
        [0] Return to main menu
        """
        print(menu)
    
    def _handle_rate_constant_calculation(self):
        """Handle calculation of rate constant from half-life."""
        print("\n===== RATE CONSTANT CALCULATOR =====")
        
        try:
            half_life = float(input("Enter half-life value: "))
            unit = input("Enter the unit of time (e.g., s, min, h): ")
            
            rate_constant = calculate_rate_constant_from_half_life(half_life)
            
            display_results_header()
            print(f"Half-life: {half_life} {unit}")
            print(f"Rate constant (k): {rate_constant:.6f} {unit}⁻¹")
            print("\nFormula used: k = ln(2) / t₁/₂")
            print(f"k = 0.693147 / {half_life} = {rate_constant:.6f} {unit}⁻¹")
        except ValueError:
            print("Error: Please enter a valid number for half-life.")
    
    def _handle_half_life_calculation(self):
        """Handle calculation of half-life from rate constant."""
        print("\n===== HALF-LIFE CALCULATOR =====")
        
        try:
            rate_constant = float(input("Enter rate constant (k) value: "))
            unit = input("Enter the unit of time (e.g., s, min, h): ")
            
            half_life = calculate_half_life_from_rate_constant(rate_constant)
            
            display_results_header()
            print(f"Rate constant (k): {rate_constant} {unit}⁻¹")
            print(f"Half-life: {half_life:.6f} {unit}")
            print("\nFormula used: t₁/₂ = ln(2) / k")
            print(f"t₁/₂ = 0.693147 / {rate_constant} = {half_life:.6f} {unit}")
        except ValueError:
            print("Error: Please enter a valid number for rate constant.")
    
    def _handle_concentration_calculation(self):
        """Handle calculation of concentration after time."""
        print("\n===== CONCENTRATION AFTER TIME CALCULATOR =====")
        
        try:
            initial_conc = float(input("Enter initial concentration [A]₀: "))
            rate_constant = float(input("Enter rate constant (k): "))
            time = float(input("Enter time elapsed: "))
            
            conc_unit = input("Enter concentration unit (e.g., M, mol/L): ")
            time_unit = input("Enter time unit (e.g., s, min, h): ")
            
            final_conc = calculate_concentration_after_time(initial_conc, rate_constant, time)
            fraction = calculate_fraction_remaining(rate_constant, time)
            
            display_results_header()
            print(f"Initial concentration [A]₀: {initial_conc} {conc_unit}")
            print(f"Rate constant (k): {rate_constant} {time_unit}⁻¹")
            print(f"Time elapsed: {time} {time_unit}")
            print(f"Final concentration [A]t: {final_conc:.6f} {conc_unit}")
            print(f"Fraction remaining: {fraction:.6f}")
            
            # Try to represent as simplified fraction if possible
            fraction_as_tuple = as_simplified_fraction(fraction)
            if fraction_as_tuple:
                num, den = fraction_as_tuple
                if den <= 100:  # Only show if denominator is reasonable
                    print(f"This fraction can be expressed as {num}/{den}")
            
            print("\nFormula used: [A]t = [A]₀ × e⁻ᵏᵗ")
            print(f"[A]t = {initial_conc} × e^(-{rate_constant} × {time})")
            print(f"[A]t = {initial_conc} × {fraction:.6f}")
            print(f"[A]t = {final_conc:.6f} {conc_unit}")
        except ValueError:
            print("Error: Please enter valid numbers for all values.")
    
    def _handle_time_calculation(self):
        """Handle calculation of time to reach a concentration."""
        print("\n===== TIME TO REACH CONCENTRATION CALCULATOR =====")
        
        try:
            initial_conc = float(input("Enter initial concentration [A]₀: "))
            final_conc = float(input("Enter final concentration [A]t: "))
            rate_constant = float(input("Enter rate constant (k): "))
            
            conc_unit = input("Enter concentration unit (e.g., M, mol/L): ")
            time_unit = input("Enter time unit (e.g., s, min, h): ")
            
            # Calculate time using ln([A]₀/[A]t) / k
            if final_conc <= 0 or final_conc >= initial_conc:
                raise ValueError("Final concentration must be positive and less than initial concentration")
            
            time = (1/rate_constant) * (initial_conc/final_conc)
            
            display_results_header()
            print(f"Initial concentration [A]₀: {initial_conc} {conc_unit}")
            print(f"Final concentration [A]t: {final_conc} {conc_unit}")
            print(f"Rate constant (k): {rate_constant} {time_unit}⁻¹")
            print(f"Time required: {time:.6f} {time_unit}")
            
            print("\nFormula used: t = ln([A]₀/[A]t) / k")
            print(f"t = ln({initial_conc}/{final_conc}) / {rate_constant}")
            print(f"t = ln({initial_conc/final_conc:.6f}) / {rate_constant}")
            print(f"t = {time:.6f} {time_unit}")
        except ValueError as e:
            print(f"Error: {e}")
    
    def _handle_complete_analysis(self):
        """Handle complete first-order kinetics analysis."""
        print("\n=== Complete First-Order Kinetics Analysis ===")
        print("This performs a complete analysis given any two parameters.")
        
        print("\nWhat information do you have? Enter values for any two of the following:")
        print("(Leave blank if unknown)")
        
        half_life_input = input("Half-life (leave blank if unknown): ")
        half_life = float(half_life_input) if half_life_input else None
        
        rate_constant_input = input("Rate constant k (leave blank if unknown): ")
        rate_constant = float(rate_constant_input) if rate_constant_input else None
        
        time_input = input("Time elapsed (leave blank if unknown): ")
        time = float(time_input) if time_input else None
        
        fraction_input = input("Fraction remaining (e.g., 0.5 for half) (leave blank if unknown): ")
        fraction_remaining = float(fraction_input) if fraction_input else None
        
        initial_conc_input = input("Initial concentration [A]₀ (leave blank if unknown): ")
        initial_concentration = float(initial_conc_input) if initial_conc_input else None
        
        final_conc_input = input("Final concentration [A]t (leave blank if unknown): ")
        final_concentration = float(final_conc_input) if final_conc_input else None
        
        unit = input("Enter the unit of time (e.g., s, min, h): ")
        conc_unit = input("Enter the concentration unit (if applicable, e.g., M, mol/L): ") if initial_concentration or final_concentration else None
        
        try:
            result = solve_first_order_kinetics(
                half_life=half_life,
                rate_constant=rate_constant,
                initial_concentration=initial_concentration,
                final_concentration=final_concentration,
                time=time,
                fraction_remaining=fraction_remaining
            )
            
            display_results_header()
            display_steps(result["steps"])
            
            print(f"\nFinal Results:")
            if result["half_life"] is not None:
                print(f"Half-life: {result['half_life']:.6f} {unit}")
            if result["rate_constant"] is not None:
                print(f"Rate constant (k): {result['rate_constant']:.6f} {unit}⁻¹")
            if result["time"] is not None:
                print(f"Time elapsed: {result['time']:.6f} {unit}")
            if result["fraction_remaining"] is not None:
                print(f"Fraction remaining: {result['fraction_remaining']:.6f}")
                if result["fraction_as_tuple"]:
                    num, den = result["fraction_as_tuple"]
                    if den <= 100:  # Only show if denominator is reasonable
                        print(f"This fraction can be expressed as {num}/{den}")
            if result["initial_concentration"] is not None:
                print(f"Initial concentration [A]₀: {result['initial_concentration']:.6f} {conc_unit}")
            if result["final_concentration"] is not None:
                print(f"Final concentration [A]t: {result['final_concentration']:.6f} {conc_unit}")
                
        except Exception as e:
            print(f"\nError: {e}")
            print("Please provide sufficient information to solve the problem.")