import xmltodict


def get_keys(x):
    dataKeys = {}

    for component in x['CDXML']['page']['stoichiometrygrid']['sgcomponent']:
        if component.get('@ComponentIsHeader', False):
            for datum in component['sgdatum']:
                dataKeys[datum['@SGPropertyType']] = datum['@SGDataValue']
    return dataKeys


def compounds(x, dataKeys,  product_ids, reagent_ids, reactant_ids):
    for component in x['CDXML']['page']['stoichiometrygrid']['sgcomponent']:
        if not component.get('@ComponentIsHeader', False):
            role = None
            if component.get('@ComponentReferenceID', "") in product_ids:
                role = 'product'
            elif component.get('@ComponentReferenceID', "") in reagent_ids:
                role = 'reagent'
            elif component.get('@ComponentReferenceID', "") in reagent_ids:
                role = 'reagent'
            dicttoyield = {
                dataKeys[datum['@SGPropertyType']]:
                datum['@SGDataValue'] for datum in component['sgdatum']
            }
            dicttoyield['role'] = role
            for k, v in dicttoyield.items():
                if "%" in k:
                    dicttoyield[k] = str(float(v.strip()) * 100)
            yield dicttoyield

 # ReactionStepReactants="5"
 # ReactionStepProducts="23"
 # ReactionStepArrows="14"
 # ReactionStepObjectsAboveArrow="17"


def parse(xml_path):
    with open(xml_path) as xfile:
        xml = xfile.read()
    x = xmltodict.parse(xml)
    if x['CDXML']['page'].get('scheme'):
        reactants = x['CDXML']['page']['scheme'][
            'step'].get('@ReactionStepReactants', None)
        reactant_ids = []
        if reactants:
            reactant_ids = reactants.split(' ')
        reagents = x['CDXML']['page']['scheme']['step'].get(
            '@ReactionStepObjectsAboveArrow', None)
        reagent_ids = []
        if reagents:
            reagent_ids = reagents.split(' ')
        products = x['CDXML']['page']['scheme'][
            'step'].get('@ReactionStepProducts', None)
        product_ids = []
        if products:
            product_ids = products.split(' ')
        keys = get_keys(x)

        data = [p for p in compounds(
            x, keys,  product_ids, reagent_ids, reactant_ids)]
        order_of_data = reactant_ids + reagent_ids + product_ids

        return {str(atomid): data[index] for index, atomid in enumerate(order_of_data)}
    print('not a reaction')
    return []
