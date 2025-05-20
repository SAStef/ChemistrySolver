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
    match = re.match(r'(\d*)\s*(.*)', term)
    coef = int(match.group(1)) if match.group(1) else 1
    formula = match.group(2).strip()
    if not formula:
        raise ValueError("Empty formula detected.")
    return coef, formula

def parse_ionic_equation(equation):
    lhs, rhs = split_ionic_equation(equation)
    reactants = [parse_term(r) for r in lhs.split('+')]
    products = [parse_term(p) for p in rhs.split('+')]
    reactants = [(c, clean_formula(f)) for c, f in reactants]
    products = [(c, clean_formula(f)) for c, f in products]
    return reactants, products

def identify_oxidation_changes(equation):
    lhs, rhs = split_ionic_equation(equation)
    reactants = [clean_formula(parse_term(r)[1]) for r in lhs.split('+')]
    products = [clean_formula(parse_term(p)[1]) for p in rhs.split('+')]
    elements = set(re.findall(r'[A-Z][a-z]?', ''.join(reactants + products)))
    changes = []
    for el in elements:
        r_ox = next((calculate_oxidation_number(f, el)['oxidation_number'] for f in reactants if el in f), None)
        p_ox = next((calculate_oxidation_number(f, el)['oxidation_number'] for f in products if el in f), None)
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
            if el['element'] in f:
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
    result.update({
        "oxidation_half": oxidation_half,
        "reduction_half": reduction_half,
        "cell_potential": cell_potential,
        "favorable": favorable
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
