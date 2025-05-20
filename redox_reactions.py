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

def clean_formula(formula):
    """
    Clean ionic formula by replacing charges with neutral notation
    
    Args:
        formula (str): Formula with possible charges like SO4^2- or Mn^2+
        
    Returns:
        str: Formula without charge notation
    """
    # Replace patterns like SO4^2- or Mn^2+
    formula = re.sub(r'(\w+)\^(\d+)([+-])', r'\1', formula)
    formula = re.sub(r'(\w+)\^([+-])(\d*)', r'\1', formula)
    formula = re.sub(r'(\w+)(\d+)([+-])', r'\1\2', formula)
    formula = re.sub(r'(\w+)([+-])(\d*)', r'\1', formula)
    return formula

def parse_ionic_equation(equation):
    """
    Parse an ionic equation, supporting charges in the format of SO4^2-, Fe^3+, etc.
    
    Args:
        equation (str): Chemical equation with ionic species
        
    Returns:
        tuple: Parsed equation suitable for balancing
    """
    # First clean the equation by removing charges
    cleaned_equation = equation
    
    # Replace charges with neutral notation for processing
    # Examples: SO4^2-, Mn^2+, H+, SO4(2-), Fe(3+)
    charge_patterns = [
        r'(\w+)\^(\d+)([+-])',   # SO4^2-
        r'(\w+)\^([+-])(\d*)',   # SO4^-
        r'(\w+)(\d+)([+-])',     # Mn2+
        r'(\w+)([+-])(\d*)'      # H+
    ]
    
    for pattern in charge_patterns:
        cleaned_equation = re.sub(pattern, lambda m: m.group(1), cleaned_equation)
    
    # Handle compounds with spaces and parentheses (like "Fe(3+)")
    cleaned_equation = re.sub(r'(\w+)\((\d+)([+-])\)', r'\1', cleaned_equation)
    cleaned_equation = re.sub(r'(\w+)\s+', r'\1 ', cleaned_equation)
    
    # Parse the cleaned equation
    try:
        return parse_chemical_equation(cleaned_equation), {}
    except Exception as e:
        # If there's an error, try a more aggressive cleaning approach
        cleaned_equation = re.sub(r'[^A-Za-z0-9\s\+\-\>\(\)]+', '', cleaned_equation)
        return parse_chemical_equation(cleaned_equation), {}

def identify_oxidation_changes(equation):
    """
    Identify elements that change oxidation states in a redox reaction.
    
    Args:
        equation (str): Balanced chemical equation
        
    Returns:
        list: Elements with their oxidation state changes
    """
    # Clean equation of any ionic charges before parsing
    clean_eq = equation
    for pattern in [r'\^(\d+)([+-])', r'\^([+-])(\d*)', r'(\d+)([+-])', r'([+-])(\d*)']:
        clean_eq = re.sub(pattern, '', clean_eq)
    
    reactants, products = parse_chemical_equation(clean_eq)
    
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
    
    # Clean the equation from ionic notation
    cleaned_equation = equation
    
    # Clean various charge notations
    charge_patterns = [
        (r'(\w+)\^(\d+)([+-])', r'\1'),  # SO4^2-
        (r'(\w+)\^([+-])(\d*)', r'\1'),  # SO4^-
        (r'(\w+)(\d+)([+-])', r'\1\2'),  # Mn2+
        (r'(\w+)([+-])(\d*)', r'\1')     # H+
    ]
    
    for pattern, replacement in charge_patterns:
        cleaned_equation = re.sub(pattern, replacement, cleaned_equation)
    
    # Step 1: Parse and initially balance the equation
    try:
        reactants, products = parse_chemical_equation(cleaned_equation)
        balanced_reactants, balanced_products = balance_equation(reactants, products)
        balanced_eq = format_balanced_equation(balanced_reactants, balanced_products)
    except Exception as e:
        # If there's an error, try a more aggressive cleaning approach
        cleaned_equation = re.sub(r'[^A-Za-z0-9\s\+\-\>\(\)]+', '', cleaned_equation)
        reactants, products = parse_chemical_equation(cleaned_equation)
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
    
    # Clean the equation and parse it
    clean_eq = equation
    for pattern in [r'\^(\d+)([+-])', r'\^([+-])(\d*)', r'(\d+)([+-])', r'([+-])(\d*)']:
        clean_eq = re.sub(pattern, '', clean_eq)
        
    reactants, products = parse_chemical_equation(clean_eq)
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
                if oxidized_element in species[1] and f"{oxidized_element}" in key:
                    oxidation_half_key = key
    
    for compound in product_compounds:
        if oxidized_element and re.search(f'({oxidized_element})', compound):
            for key in STANDARD_REDUCTION_POTENTIALS:
                species = key.split('/')
                if oxidized_element in species[0] and f"{oxidized_element}" in key:
                    oxidation_half_key = key
    
    # Look for reduced element in reactants and products
    for compound in reactant_compounds:
        if reduced_element and re.search(f'({reduced_element})', compound):
            for key in STANDARD_REDUCTION_POTENTIALS:
                species = key.split('/')
                if reduced_element in species[0] and f"{reduced_element}" in key:
                    reduction_half_key = key
    
    for compound in product_compounds:
        if reduced_element and re.search(f'({reduced_element})', compound):
            for key in STANDARD_REDUCTION_POTENTIALS:
                species = key.split('/')
                if reduced_element in species[1] and f"{reduced_element}" in key:
                    reduction_half_key = key
    
    # Direct lookup for common half-reactions
    reaction_lookup = {
        "Zn + Cd2+ -> Zn2+ + Cd": {
            "oxidation": "Zn2+/Zn",
            "reduction": "Cd2+/Cd"
        },
        "Cu + 2Ag+ -> Cu2+ + 2Ag": {
            "oxidation": "Cu2+/Cu",
            "reduction": "Ag+/Ag"
        },
        # Add more common reactions as needed
    }
    
    # Try direct lookup first
    clean_equation = re.sub(r'\s+', ' ', equation).strip()
    for key in reaction_lookup:
        # Create a pattern-matching version that ignores whitespace and arrow style
        pattern_key = re.sub(r'\s+', '\\s*', key)
        pattern_key = pattern_key.replace("->", "(?:->|→|-->)")
        pattern_key = pattern_key.replace("+", "\\+")
        
        if re.search(pattern_key, clean_equation, re.IGNORECASE):
            oxidation_half_key = reaction_lookup[key]["oxidation"]
            reduction_half_key = reaction_lookup[key]["reduction"]
            break
    
    # Alternative approach: Parse the reaction and identify the redox pairs
    if not oxidation_half_key or not reduction_half_key:
        for half_key in STANDARD_REDUCTION_POTENTIALS:
            species = half_key.split('/')
            
            # Clean species names for matching
            clean_species0 = re.sub(r'(\d+)([+-])', r'\1', species[0])
            clean_species1 = species[1]
            
            # Check if both species from a half-reaction are in our equation
            if (clean_species0 in clean_eq and clean_species1 in clean_eq):
                # Determine if this is oxidation or reduction in our reaction
                if clean_species0 in ''.join(product_compounds) and clean_species1 in ''.join(reactant_compounds):
                    oxidation_half_key = half_key
                elif clean_species0 in ''.join(reactant_compounds) and clean_species1 in ''.join(product_compounds):
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
        # Clean the equation of ionic charges
        clean_eq = equation
        for pattern in [r'\^(\d+)([+-])', r'\^([+-])(\d*)', r'(\d+)([+-])', r'([+-])(\d*)']:
            clean_eq = re.sub(pattern, '', clean_eq)
            
        # Clean the compounds as well
        clean_compound1 = clean_formula(compound1)
        clean_compound2 = clean_formula(compound2)
        
        reactants, products = parse_chemical_equation(clean_eq)
        balanced_reactants, balanced_products = balance_equation(reactants, products)
        
        # Combine all compounds from both sides
        all_compounds = {formula: coef for coef, formula in balanced_reactants + balanced_products}
        
        # Check if both compounds exist in the equation
        if clean_compound1 in all_compounds and clean_compound2 in all_compounds:
            ratio = all_compounds[clean_compound1] / all_compounds[clean_compound2]
            return ratio
        
        return None
    except Exception as e:
        print(f"Error in finding molar ratio: {str(e)}")
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
    
    # Clean equation to standardize compound names
    clean_eq = result["balanced_equation"]
    
    # Check if SO2 is involved
    if "SO2" in clean_eq:
        # Find molar ratio between SO2 and H+
        if "H" in clean_eq:  # We'll use H as a simplification of H+
            ratio = find_molar_ratio(clean_eq, "SO2", "H")
            if ratio:
                molar_ratio = (1, 1/ratio)  # Convert to SO2:H+ ratio
    
    # Add the molar ratio to the result
    result["molar_ratio"] = molar_ratio
    
    return result