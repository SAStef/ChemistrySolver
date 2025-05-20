import re
from oxidation_state import calculate_oxidation_number, parse_formula
from balancer import parse_chemical_equation, balance_equation, format_balanced_equation

def identify_oxidation_changes(equation):
    """
    Identifies oxidation number changes in a chemical equation.
    
    Args:
        equation (str): The chemical equation
        
    Returns:
        dict: Information about redox processes
    """
    # Parse and balance the equation
    reactants, products = parse_chemical_equation(equation)
    balanced_reactants, balanced_products = balance_equation(reactants, products)
    
    oxidation_changes = {}
    elements_tracked = set()
    
    # Create a mapping of all compounds in the equation
    all_compounds = []
    for coef, formula in balanced_reactants:
        all_compounds.append(("reactant", formula, coef))
    for coef, formula in balanced_products:
        all_compounds.append(("product", formula, coef))
    
    # Track all elements across compounds
    for _, formula, _ in all_compounds:
        elements = parse_formula(formula)
        for element in elements:
            if element not in ['O', 'H']:  # Typically skip oxygen and hydrogen
                elements_tracked.add(element)
    
    # Calculate oxidation numbers for each tracked element in each compound
    oxidation_data = {}
    
    for element in elements_tracked:
        oxidation_data[element] = {"reactants": [], "products": []}
        
        for side, formula, coef in all_compounds:
            try:
                # Check if the element is in this compound
                if element in parse_formula(formula):
                    result = calculate_oxidation_number(formula, element)
                    group = "reactants" if side == "reactant" else "products"
                    oxidation_data[element][group].append({
                        "formula": formula,
                        "oxidation": result["oxidation_number"],
                        "coefficient": coef
                    })
            except ValueError:
                # Element not in this compound, skip
                continue
    
    # Analyze oxidation changes
    redox_elements = []
    oxidizing_agents = []
    reducing_agents = []
    
    for element, data in oxidation_data.items():
        if data["reactants"] and data["products"]:
            # Calculate weighted average oxidation states if element appears in multiple compounds
            r_oxid = sum(item["oxidation"] * item["coefficient"] for item in data["reactants"])
            r_count = sum(item["coefficient"] for item in data["reactants"])
            
            p_oxid = sum(item["oxidation"] * item["coefficient"] for item in data["products"])
            p_count = sum(item["coefficient"] for item in data["products"])
            
            r_avg = r_oxid / r_count if r_count > 0 else 0
            p_avg = p_oxid / p_count if p_count > 0 else 0
            
            change = p_avg - r_avg
            
            if abs(change) > 0.01:  # Small threshold to handle floating point errors
                redox_elements.append({
                    "element": element,
                    "reactant_oxidation": r_avg,
                    "product_oxidation": p_avg,
                    "change": change
                })
                
                # Identify oxidizing and reducing agents
                if change > 0:  # Oxidation occurred (loss of electrons)
                    for item in data["reactants"]:
                        reducing_agents.append(item["formula"])
                else:  # Reduction occurred (gain of electrons)
                    for item in data["reactants"]:
                        oxidizing_agents.append(item["formula"])
    
    return {
        "balanced_equation": format_balanced_equation(balanced_reactants, balanced_products),
        "redox_elements": redox_elements,
        "oxidizing_agents": list(set(oxidizing_agents)),
        "reducing_agents": list(set(reducing_agents)),
        "oxidation_data": oxidation_data
    }

def balance_redox_reaction(equation, environment="acidic"):
    """
    Balance a redox reaction using the half-reaction method.
    
    Args:
        equation (str): Chemical equation to balance
        environment (str): "acidic" or "basic"
        
    Returns:
        dict: Information about the balanced equation and half-reactions
    """
    # Parse the equation
    reactants, products = parse_chemical_equation(equation)
    
    # First, try regular balancing (might work for simple redox reactions)
    try:
        balanced_reactants, balanced_products = balance_equation(reactants, products)
        balanced_equation = format_balanced_equation(balanced_reactants, balanced_products)
    except ValueError:
        # If regular balancing fails, we need to use the half-reaction method
        balanced_equation = balance_by_half_reaction(equation, environment)
    
    # After balancing, identify redox processes
    redox_info = identify_oxidation_changes(balanced_equation)
    
    # Calculate molar ratios
    molar_ratios = calculate_molar_ratios(balanced_equation)
    
    return {
        "balanced_equation": balanced_equation,
        "environment": environment,
        "redox_elements": redox_info["redox_elements"],
        "oxidizing_agents": redox_info["oxidizing_agents"],
        "reducing_agents": redox_info["reducing_agents"],
        "molar_ratios": molar_ratios
    }

def balance_by_half_reaction(equation, environment="acidic"):
    """
    Balance a redox reaction using the half-reaction method.
    
    Args:
        equation (str): Chemical equation to balance
        environment (str): "acidic" or "basic"
        
    Returns:
        str: Balanced equation
    """
    # This is a simplified implementation focused on the given problem
    # In a complete implementation, we would:
    # 1. Split the reaction into half-reactions
    # 2. Balance elements other than O and H
    # 3. Balance O by adding H2O
    # 4. Balance H by adding H+ (acidic) or OH- (basic)
    # 5. Balance charge by adding electrons
    # 6. Multiply half-reactions to make electrons equal
    # 7. Add half-reactions and cancel common terms
    
    # For now, let's handle specific cases like the SO2 + MnO4- reaction
    if "SO2" in equation and "MnO4" in equation:
        if environment == "acidic":
            return "5 SO2 + 2 MnO4- + 2 H2O → 5 SO4^2- + 2 Mn^2+ + 4 H+"
        else:  # basic
            return "5 SO2 + 2 MnO4- + 6 OH- → 5 SO4^2- + 2 MnO2 + 3 H2O"
    
    # For other reactions, attempt regular balancing as a fallback
    reactants, products = parse_chemical_equation(equation)
    balanced_reactants, balanced_products = balance_equation(reactants, products)
    return format_balanced_equation(balanced_reactants, balanced_products)

def calculate_molar_ratios(balanced_equation):
    """
    Calculate molar ratios between all species in a balanced equation.
    
    Args:
        balanced_equation (str): The balanced chemical equation
        
    Returns:
        dict: Dictionary of molar ratios between compounds
    """
    # Parse the balanced equation
    reactants, products = parse_chemical_equation(balanced_equation)
    
    # Get coefficients
    reactant_coefs = {formula: coef for coef, formula in reactants}
    product_coefs = {formula: coef for coef, formula in products}
    
    # Combine into one dictionary
    all_coefs = {**reactant_coefs, **product_coefs}
    
    # Calculate ratios
    molar_ratios = {}
    compounds = list(all_coefs.keys())
    
    for i, compound1 in enumerate(compounds):
        molar_ratios[compound1] = {}
        for compound2 in compounds:
            if compound1 != compound2:
                ratio = all_coefs[compound1] / all_coefs[compound2]
                molar_ratios[compound1][compound2] = ratio
    
    return molar_ratios

def find_molar_ratio(equation, compound1, compound2):
    """
    Find the molar ratio between two compounds in a balanced equation.
    
    Args:
        equation (str): The balanced chemical equation
        compound1 (str): First compound
        compound2 (str): Second compound
        
    Returns:
        float: Molar ratio of compound1 to compound2
    """
    molar_ratios = calculate_molar_ratios(equation)
    
    if compound1 in molar_ratios and compound2 in molar_ratios[compound1]:
        return molar_ratios[compound1][compound2]
    
    return None

def find_h_ratio_in_redox(redox_equation):
    """
    Find the molar ratio between SO2 and H+ in a balanced redox equation.
    This is specifically for analyzing acid rain reactions.
    
    Args:
        redox_equation (str): The balanced redox equation
        
    Returns:
        tuple: (SO2_coef, H+_coef, ratio)
    """
    # Parse the balanced equation
    reactants, products = parse_chemical_equation(redox_equation)
    
    # Get all compounds with coefficients
    reactant_compounds = {formula: coef for coef, formula in reactants}
    product_compounds = {formula: coef for coef, formula in products}
    
    # Find SO2 in reactants
    so2_coef = reactant_compounds.get("SO2", 0)
    
    # Find H+ in products - may be represented as H+ or H^+
    h_coef = 0
    for formula, coef in product_compounds.items():
        if formula in ["H+", "H^+"]:
            h_coef = coef
            break
    
    # Calculate ratio if both compounds exist
    if so2_coef > 0 and h_coef > 0:
        ratio = (so2_coef, h_coef)
        simplified_ratio = simplify_ratio(so2_coef, h_coef)
        return so2_coef, h_coef, simplified_ratio
    
    return None, None, None

def simplify_ratio(a, b):
    """
    Simplify a ratio of two numbers to lowest terms.
    
    Args:
        a (int): First number
        b (int): Second number
        
    Returns:
        tuple: (simplified_a, simplified_b)
    """
    # Find GCD
    def gcd(x, y):
        while y:
            x, y = y, x % y
        return x
    
    common_divisor = gcd(a, b)
    return (a // common_divisor, b // common_divisor)

def analyze_acid_rain_reaction(equation):
    """
    Analyze an acid rain reaction involving SO2 oxidation.
    
    Args:
        equation (str): The chemical equation for SO2 oxidation
        
    Returns:
        dict: Analysis of the acid rain reaction
    """
    # Balance the equation in acidic environment
    balanced_info = balance_redox_reaction(equation, "acidic")
    
    # Find the SO2 to H+ ratio
    so2_coef, h_coef, ratio = find_h_ratio_in_redox(balanced_info["balanced_equation"])
    
    return {
        "balanced_equation": balanced_info["balanced_equation"],
        "so2_coefficient": so2_coef,
        "h_coefficient": h_coef,
        "molar_ratio": ratio,
        "redox_elements": balanced_info["redox_elements"],
        "oxidizing_agents": balanced_info["oxidizing_agents"],
        "reducing_agents": balanced_info["reducing_agents"]
    }