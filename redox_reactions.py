import re
from oxidation_state import calculate_oxidation_number

STANDARD_REDUCTION_POTENTIALS = {
    "Zn2+/Zn": -0.76,
    "Cd2+/Cd": -0.40,
    "Cu2+/Cu": 0.34,
    "Ag+/Ag": 0.80,
    # Extend this dictionary with more half-reactions as needed
}

def clean_formula(formula):
    return re.sub(r'([A-Za-z0-9]+)(\^\d+[+-]|\^?[+-]\d*|\(\d+[+-]\)|[+-])', r'\1', formula)

def split_ionic_equation(equation):
    # Only split on arrows — do not remove charges at this stage
    for sep in ['→', '->', '-->']:
        if sep in equation:
            idx = equation.find(sep)
            return [equation[:idx].strip(), equation[idx + len(sep):].strip()]
    raise ValueError("Invalid equation format. Use '->' or '→'.")

def parse_term(term):
    # Fix the issue with handling coefficient and formula
    term = term.strip()
    match = re.match(r'^(\d*)\s*(.+)$', term)
    
    if not match:
        raise ValueError(f"Invalid term format: {term}")
    
    coef = int(match.group(1)) if match.group(1) else 1
    formula = match.group(2).strip()
    
    if not formula:
        raise ValueError(f"Empty formula detected in term: '{term}'")
        
    return coef, formula

def parse_ionic_equation(equation):
    lhs, rhs = split_ionic_equation(equation)
    
    # Handle case of no spaces around plus sign
    lhs = re.sub(r'(\S)\+(\S)', r'\1 + \2', lhs)
    rhs = re.sub(r'(\S)\+(\S)', r'\1 + \2', rhs)
    
    reactants = [parse_term(r.strip()) for r in lhs.split('+')]
    products = [parse_term(p.strip()) for p in rhs.split('+')]
    
    reactants = [(c, clean_formula(f)) for c, f in reactants]
    products = [(c, clean_formula(f)) for c, f in products]
    
    return reactants, products

def identify_oxidation_changes(equation):
    reactants, products = parse_ionic_equation(equation)
    reactant_formulas = [f for _, f in reactants]
    product_formulas = [f for _, f in products]
    
    elements = set()
    for formula in reactant_formulas + product_formulas:
        elements.update(re.findall(r'[A-Z][a-z]?', formula))
    
    changes = []
    for el in elements:
        r_ox = None
        p_ox = None
        
        # Find oxidation state in reactants
        for formula in reactant_formulas:
            if re.search(f'[^A-Za-z]{el}|^{el}', formula):
                try:
                    r_ox = calculate_oxidation_number(formula, el)['oxidation_number']
                    break
                except:
                    pass
        
        # Find oxidation state in products
        for formula in product_formulas:
            if re.search(f'[^A-Za-z]{el}|^{el}', formula):
                try:
                    p_ox = calculate_oxidation_number(formula, el)['oxidation_number']
                    break
                except:
                    pass
        
        if r_ox is not None and p_ox is not None and r_ox != p_ox:
            changes.append({
                "element": el,
                "reactant_oxidation": r_ox,
                "product_oxidation": p_ox,
                "change": p_ox - r_ox
            })
    
    return changes

def balance_redox_reaction(equation, environment="acidic"):
    reactants, products = parse_ionic_equation(equation)
    redox_elements = identify_oxidation_changes(equation)
    oxidizing_agents = []
    reducing_agents = []
    
    for el in redox_elements:
        for _, f in reactants:
            if re.search(f'[^A-Za-z]{el["element"]}|^{el["element"]}', f):
                if el['change'] > 0:
                    reducing_agents.append(f)
                elif el['change'] < 0:
                    oxidizing_agents.append(f)
    
    # Simple formatted string for visual output
    balanced_eq = ' + '.join(f"{c if c > 1 else ''}{f}" for c, f in reactants)
    balanced_eq += ' -> '
    balanced_eq += ' + '.join(f"{c if c > 1 else ''}{f}" for c, f in products)
    
    return {
        "balanced_equation": balanced_eq,
        "environment": environment,
        "redox_elements": redox_elements,
        "oxidizing_agents": list(set(oxidizing_agents)),
        "reducing_agents": list(set(reducing_agents))
    }

def determine_reaction_favorability(equation):
    result = balance_redox_reaction(equation)
    redox_elements = result['redox_elements']
    
    oxidized = next((el for el in redox_elements if el['change'] > 0), None)
    reduced = next((el for el in redox_elements if el['change'] < 0), None)
    
    oxidation_half, reduction_half = None, None
    
    # Find matching half-reactions in our database
    for key in STANDARD_REDUCTION_POTENTIALS:
        sp = key.split('/')
        if oxidized and oxidized['element'] in sp:
            oxidation_half = key
        if reduced and reduced['element'] in sp:
            reduction_half = key
    
    ox_pot = STANDARD_REDUCTION_POTENTIALS.get(oxidation_half)
    red_pot = STANDARD_REDUCTION_POTENTIALS.get(reduction_half)
    
    cell_potential = None
    if ox_pot is not None and red_pot is not None:
        cell_potential = red_pot - ox_pot
    
    favorable = cell_potential is not None and cell_potential > 0
    
    explanation = ""
    if cell_potential is not None:
        if favorable:
            explanation = f"The reaction is favorable because the cell potential (E°cell) is positive ({cell_potential:.2f} V). "
            explanation += f"This indicates that the reaction will proceed spontaneously under standard conditions."
        else:
            explanation += f"The reaction is not favorable because the cell potential (E°cell) is negative ({cell_potential:.2f} V). "
            explanation += f"The reverse reaction would be spontaneous under standard conditions."
    
    result.update({
        "oxidation_half": oxidation_half,
        "reduction_half": reduction_half,
        "cell_potential": cell_potential,
        "favorable": favorable,
        "explanation": explanation
    })
    
    return result

def find_molar_ratio(equation, compound1, compound2):
    try:
        reactants, products = parse_ionic_equation(equation)
        all_compounds = {f: c for c, f in reactants + products}
        
        c1 = clean_formula(compound1)
        c2 = clean_formula(compound2)
        
        if c1 in all_compounds and c2 in all_compounds:
            return all_compounds[c1] / all_compounds[c2]
        
        return None
    except Exception as e:
        print(f"Error in finding molar ratio: {str(e)}")
        return None