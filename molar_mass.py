try:
    import molmass
except ImportError:
    molmass = None

def check_molmass():
    if molmass is None:
        raise ImportError("molmass is required. Install it with `pip install molmass`.")

def calculate_molar_mass(formula):
    try:
        check_molmass()
        f = molmass.Formula(formula)
        comp = [{
            'element': s,
            'count': i.count,
            'atomic_mass': i.mass / i.count,
            'contribution': i.mass
        } for s, i in f.composition().items()]

        return {
            'formula': str(f),
            'hill_notation': f.formula,
            'molar_mass': f.mass,
            'monoisotopic_mass': getattr(f.isotope, 'mass', f.mass),
            'composition': comp,
            'success': True
        }
    except Exception as e:
        return {'success': False, 'error': str(e)}
