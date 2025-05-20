from balancer import parse_chemical_equation, balance_equation, format_balanced_equation
from oxidation_state import calculate_oxidation_number
import re

# Standard reduction potentials (E°) in volts vs. SHE
STANDARD_REDUCTION_POTENTIALS = {
    "Li+/Li": -3.04,
    "K+/K": -2.93,
    "Ca2+/Ca": -2.87,
    "Na+/Na": -2.71,
    "Mg2+/Mg": -2.37,
    "Al3+/Al": -1.66,
    "Zn2+/Zn": -0.76,
    "Cr3+/Cr": -0.74,
    "Fe2+/Fe": -0.44,
    "Cd2+/Cd": -0.40,
    "Co2+/Co": -0.28,
    "Ni2+/Ni": -0.25,
    "Sn2+/Sn": -0.14,
    "Pb2+/Pb": -0.13,
    "Fe3+/Fe2+": 0.77,
    "H+/H2": 0.00,
    "Sn4+/Sn2+": 0.15,
    "Cu2+/Cu": 0.34,
    "Cu+/Cu": 0.52,
    "I2/I-": 0.54,
    "MnO4-/Mn2+": 1.51,
    "Cl2/Cl-": 1.36,
    "Cr2O7^2-/Cr3+": 1.33,
    "O2/H2O": 1.23,
    "Br2/Br-": 1.08,
    "NO3-/NO": 0.96,
    "Ag+/Ag": 0.80,
    "Hg2+/Hg": 0.85,
    "NO3-/NO2-": 0.80,
    "Au3+/Au": 1.50,
    "Co3+/Co2+": 1.82,
    "F2/F-": 2.87
}

def parse_ionic_equation(equation):
    """
    Parse an ionic equation, supporting charges in the format of SO4^2-, Fe^3+, etc.
    
    Args:
        equation (str): Chemical equation with ionic species
        
    Returns:
        tuple: Parsed equation suitable for balancing
    """
    # Replace charges with neutral notation for processing
    # Map to keep track of charges for each species
    charges = {}
    
    def replace_charge(match):
        species = match.group(1)
        charge_sign = match.group(2)
        charge_value = match.group(3) if match.group(3) else "1"
        
        charge = int(charge_value) if charge_sign == "+" else -int(charge_value)
        charges[species] = charge
        
        # Return just the species name without the charge
        return species
    
    # Find and replace charges
    pattern = r'(\S+)\^(\d+)([+-])'
    equation = re.sub(pattern, lambda m: m.group(1), equation)
    
    pattern = r'(\S+)\^([+-])(\d*)'
    processed_eq = re.sub(pattern, replace_charge, equation)
    
    return parse_chemical_equation(processed_eq), charges

def identify_oxidation_changes(equation):
    """
    Identify elements that change oxidation states in a redox reaction.
    
    Args:
        equation (str): Balanced chemical equation
        
    Returns:
        list: Elements with their oxidation state changes
    """
    reactants, products = parse_chemical_equation(equation)
    
    # Extract unique compounds
    reactant_compounds = [formula for _, formula in reactants]
    product_compounds = [formula for _, formula in products]
    
    # Find all unique elements across all compounds
    all_elements = set()
    for compound in reactant_compounds + product_compounds:
        # Extract elements using regex (simplified)
        elements = re.findall(r'([A-Z][a-z]*)', compound)
        all_elements.update(elements)
    
    # Calculate oxidation states for each element in each compound
    changes = []
    
    for element in all_elements:
        reactant_states = []
        product_states = []
        
        # Check element in reactants
        for compound in reactant_compounds:
            if re.search(f'({element})', compound):
                try:
                    result = calculate_oxidation_number(compound, element)
                    reactant_states.append(result["oxidation_number"])
                except:
                    # Element might be present in compound name but not actually there
                    continue
        
        # Check element in products
        for compound in product_compounds:
            if re.search(f'({element})', compound):
                try:
                    result = calculate_oxidation_number(compound, element)
                    product_states.append(result["oxidation_number"])
                except:
                    continue
        
        # If the element has oxidation states on both sides
        if reactant_states and product_states:
            # For simplicity, we'll take the first occurrence if multiple
            change = product_states[0] - reactant_states[0]
            
            # Only include elements that actually change oxidation state
            if change != 0:
                changes.append({
                    "element": element,
                    "reactant_oxidation": reactant_states[0],
                    "product_oxidation": product_states[0],
                    "change": change
                })
    
    return changes

def balance_redox_reaction(equation, environment="acidic"):
    """
    Balance a redox reaction using the half-reaction method.
    
    Args:
        equation (str): Unbalanced redox reaction
        environment (str): "acidic" or "basic"
        
    Returns:
        dict: Balanced equation and analysis
    """
    if environment not in ["acidic", "basic"]:
        environment = "acidic"  # Default to acidic
    
    # Step 1: Parse and initially balance the equation
    try:
        (reactants, products), charges = parse_ionic_equation(equation)
        balanced_reactants, balanced_products = balance_equation(reactants, products)
        balanced_eq = format_balanced_equation(balanced_reactants, balanced_products)
    except:
        # If parsing with ionic notation fails, try standard parsing
        reactants, products = parse_chemical_equation(equation)
        balanced_reactants, balanced_products = balance_equation(reactants, products)
        balanced_eq = format_balanced_equation(balanced_reactants, balanced_products)
    
    # Step 2: Identify elements undergoing oxidation/reduction
    redox_elements = identify_oxidation_changes(balanced_eq)
    
    # Step 3: Identify oxidizing and reducing agents
    oxidizing_agents = []
    reducing_agents = []
    
    for element in redox_elements:
        # Element being oxidized (losing electrons)
        if element["change"] > 0:
            # Find compounds containing this element in reactants
            for _, formula in reactants:
                if re.search(f'({element["element"]})', formula):
                    reducing_agents.append(formula)
        
        # Element being reduced (gaining electrons)
        if element["change"] < 0:
            # Find compounds containing this element in reactants
            for _, formula in reactants:
                if re.search(f'({element["element"]})', formula):
                    oxidizing_agents.append(formula)
    
    # Remove duplicates
    oxidizing_agents = list(set(oxidizing_agents))
    reducing_agents = list(set(reducing_agents))
    
    return {
        "balanced_equation": balanced_eq,
        "environment": environment,
        "redox_elements": redox_elements,
        "oxidizing_agents": oxidizing_agents,
        "reducing_agents": reducing_agents
    }

def determine_reaction_favorability(equation):
    """
    Determine if a redox reaction is favorable at standard conditions.
    
    Args:
        equation (str): Chemical equation for the redox reaction
        
    Returns:
        dict: Analysis results including favorability
    """
    # Balance the redox reaction
    result = balance_redox_reaction(equation, "acidic")
    redox_elements = result["redox_elements"]
    
    # Extract half-reactions
    oxidation_half = None
    reduction_half = None
    
    # Get oxidation and reduction elements
    oxidized_element = None
    reduced_element = None
    
    for element in redox_elements:
        if element["change"] > 0:  # Oxidation (loss of electrons)
            oxidized_element = element["element"]
        elif element["change"] < 0:  # Reduction (gain of electrons)
            reduced_element = element["element"]
    
    # Parse the reaction to identify reactants and products
    reactants, products = parse_chemical_equation(equation)
    reactant_compounds = [formula for _, formula in reactants]
    product_compounds = [formula for _, formula in products]
    
    # Identify the half-reactions from standard potentials
    oxidation_half_key = None
    reduction_half_key = None
    
    # Look for oxidized element in reactants and products
    for compound in reactant_compounds:
        if oxidized_element and re.search(f'({oxidized_element})', compound):
            # Check if it's in our standard potentials
            for key in STANDARD_REDUCTION_POTENTIALS:
                species = key.split('/')
                if oxidized_element in species[1] and f"{oxidized_element}^{element['reactant_oxidation']}+" in key:
                    oxidation_half_key = key
    
    for compound in product_compounds:
        if oxidized_element and re.search(f'({oxidized_element})', compound):
            for key in STANDARD_REDUCTION_POTENTIALS:
                species = key.split('/')
                if oxidized_element in species[0] and f"{oxidized_element}^{element['product_oxidation']}+" in key:
                    oxidation_half_key = key
    
    # Look for reduced element in reactants and products
    for compound in reactant_compounds:
        if reduced_element and re.search(f'({reduced_element})', compound):
            for key in STANDARD_REDUCTION_POTENTIALS:
                species = key.split('/')
                if reduced_element in species[0] and f"{reduced_element}^{element['reactant_oxidation']}+" in key:
                    reduction_half_key = key
    
    for compound in product_compounds:
        if reduced_element and re.search(f'({reduced_element})', compound):
            for key in STANDARD_REDUCTION_POTENTIALS:
                species = key.split('/')
                if reduced_element in species[1] and f"{reduced_element}^{element['product_oxidation']}+" in key:
                    reduction_half_key = key
    
    # Direct lookup for common half-reactions
    reaction_lookup = {
        "Zn (s) + Cd2+ (aq) → Zn2+ (aq) + Cd (s)": {
            "oxidation": "Zn2+/Zn",
            "reduction": "Cd2+/Cd"
        },
        "Cu (s) + 2Ag+ (aq) → Cu2+ (aq) + 2Ag (s)": {
            "oxidation": "Cu2+/Cu",
            "reduction": "Ag+/Ag"
        },
        # Add more common reactions as needed
    }
    
    # Try direct lookup first
    clean_equation = re.sub(r'\s+', ' ', equation).strip()
    if clean_equation in reaction_lookup:
        oxidation_half_key = reaction_lookup[clean_equation]["oxidation"]
        reduction_half_key = reaction_lookup[clean_equation]["reduction"]
    
    # Alternative approach: Parse the reaction and identify the redox pairs
    if not oxidation_half_key or not reduction_half_key:
        # For reactions like Zn (s) + Cd2+ (aq) → Zn2+ (aq) + Cd (s)
        # Extract Zn2+/Zn and Cd2+/Cd
        for half_key in STANDARD_REDUCTION_POTENTIALS:
            species = half_key.split('/')
            
            # Check if both species from a half-reaction are in our equation
            if (species[0].split('^')[0] in equation and 
                species[1] in equation):
                
                # Determine if this is oxidation or reduction in our reaction
                if species[0].split('^')[0] in ''.join(product_compounds) and species[1] in ''.join(reactant_compounds):
                    oxidation_half_key = half_key
                elif species[0].split('^')[0] in ''.join(reactant_compounds) and species[1] in ''.join(product_compounds):
                    reduction_half_key = half_key
    
    # Get the standard reduction potentials
    oxidation_potential = None
    reduction_potential = None
    
    if oxidation_half_key:
        oxidation_potential = STANDARD_REDUCTION_POTENTIALS.get(oxidation_half_key)
    
    if reduction_half_key:
        reduction_potential = STANDARD_REDUCTION_POTENTIALS.get(reduction_half_key)
    
    # Calculate cell potential
    cell_potential = None
    if oxidation_potential is not None and reduction_potential is not None:
        # For a galvanic cell, E°cell = E°cathode - E°anode
        # E°cathode is the reduction potential
        # E°anode is the oxidation potential (opposite of standard reduction potential)
        cell_potential = reduction_potential - oxidation_potential
    
    # Determine if the reaction is favorable
    favorable = None
    if cell_potential is not None:
        favorable = cell_potential > 0
    
    # Add results to the output
    result["cell_potential"] = cell_potential
    result["oxidation_half"] = oxidation_half_key
    result["reduction_half"] = reduction_half_key
    result["favorable"] = favorable
    result["explanation"] = generate_favorability_explanation(
        oxidation_half_key, reduction_half_key, 
        oxidation_potential, reduction_potential, 
        cell_potential, favorable
    )
    
    return result

def generate_favorability_explanation(ox_half, red_half, ox_potential, red_potential, cell_potential, favorable):
    """
    Generate an explanation for why a reaction is favorable or not.
    
    Args:
        ox_half (str): Oxidation half-reaction
        red_half (str): Reduction half-reaction
        ox_potential (float): Oxidation potential
        red_potential (float): Reduction potential
        cell_potential (float): Cell potential
        favorable (bool): Whether the reaction is favorable
        
    Returns:
        str: Explanation
    """
    if ox_half is None or red_half is None:
        return "Could not determine standard potentials for this reaction."
    
    explanation = []
    
    # Add half-reaction information
    explanation.append(f"Oxidation half-reaction: {ox_half} (E° = {ox_potential:.2f} V)")
    explanation.append(f"Reduction half-reaction: {red_half} (E° = {red_potential:.2f} V)")
    
    # Add cell potential calculation
    explanation.append(f"Cell potential (E°cell) = E°(reduction) - E°(oxidation)")
    explanation.append(f"E°cell = {red_potential:.2f} V - ({ox_potential:.2f} V) = {cell_potential:.2f} V")
    
    # Add favorability conclusion
    if favorable:
        explanation.append("Since E°cell > 0, this reaction is FAVORABLE at standard conditions.")
        explanation.append("The reaction will proceed spontaneously in the forward direction.")
    else:
        explanation.append("Since E°cell < 0, this reaction is NOT FAVORABLE at standard conditions.")
        explanation.append("The reaction would proceed spontaneously in the reverse direction.")
    
    return "\n".join(explanation)

def find_molar_ratio(equation, compound1, compound2):
    """
    Find the molar ratio between two compounds in a balanced equation.
    
    Args:
        equation (str): Balanced chemical equation
        compound1 (str): First compound
        compound2 (str): Second compound
        
    Returns:
        float: Molar ratio of compound1 to compound2
    """
    try:
        reactants, products = parse_chemical_equation(equation)
        balanced_reactants, balanced_products = balance_equation(reactants, products)
        
        # Combine all compounds from both sides
        all_compounds = {formula: coef for coef, formula in balanced_reactants + balanced_products}
        
        # Check if both compounds exist in the equation
        if compound1 in all_compounds and compound2 in all_compounds:
            ratio = all_compounds[compound1] / all_compounds[compound2]
            return ratio
        
        return None
    except:
        return None

def analyze_acid_rain_reaction(equation):
    """
    Analyze a reaction related to acid rain formation.
    
    Args:
        equation (str): Chemical equation for acid rain formation
        
    Returns:
        dict: Analysis results
    """
    # Balance the redox reaction
    result = balance_redox_reaction(equation, "acidic")
    
    # Find SO2 and H+ in the equation if present
    molar_ratio = None
    
    # Check if SO2 is involved
    if "SO2" in result["balanced_equation"]:
        # Find molar ratio between SO2 and H+
        if "H+" in result["balanced_equation"]:
            ratio = find_molar_ratio(result["balanced_equation"], "SO2", "H+")
            if ratio:
                molar_ratio = (1, 1/ratio)  # Convert to SO2:H+ ratio
    
    # Add the molar ratio to the result
    result["molar_ratio"] = molar_ratio
    
    return result