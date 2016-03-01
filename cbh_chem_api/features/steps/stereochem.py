from behave import given, when, then
import json
import os
from django.contrib.auth.models import User, Group
from pybel import readfile


@given('I have valid stereochem compound with one stereocentre')
def step(context):
    context.inchi = "InChI=1S/C10H12O2/c1-3-8-4-7(2)5-9(6-8)10(11)12/h4-6H,3H2,1-2H3,(H,11,12)"
    context.post_data["ctab"] = """


 12 12  0  0  0  0            999 V2000
    0.4330    0.2500    0.0000 C   0  0  0  0  0  0
    1.2990   -0.2500    0.0000 C   0  0  0  0  0  0
    1.2990   -1.2500    0.0000 C   0  0  0  0  0  0
    0.4330   -1.7500    0.0000 C   0  0  0  0  0  0
   -0.4330   -1.2500    0.0000 C   0  0  0  0  0  0
   -0.4330   -0.2500    0.0000 C   0  0  0  0  0  0
    0.4330    1.2500    0.0000 C   0  0  0  0  0  0
    1.2990    1.7500    0.0000 O   0  0  0  0  0  0
   -0.4330    1.7500    0.0000 O   0  0  0  0  0  0
   -1.2990   -1.7500    0.0000 C   0  0  0  0  0  0
   -2.1651   -1.2500    0.0000 C   0  0  0  0  0  0
    2.1651   -1.7500    0.0000 C   0  0  0  0  0  0
  1  21  0     0  0
  2  32  0     0  0
  3  41  0     0  0
  4  52  0     0  0
  5  61  0     0  0
  6  12  0     0  0
  1  71  0     0  0
  7  81  0     0  0
  7  92  0     0  0
  5 101  0     0  0
 10 111  1     0  0
  3 121  0     0  0
M  END"""


@given('I have valid stereochem compound with multiple stereocentres')
def step(context):
    context.inchi = "InChI=1S/C11H14O2/c1-3-8-5-9(4-2)7-10(6-8)11(12)13/h5-7H,3-4H2,1-2H3,(H,12,13)"
    context.post_data["ctab"] = """


 13 13  0  0  0  0            999 V2000
   -0.0000    0.2500    0.0000 C   0  0  0  0  0  0
    0.8660   -0.2500    0.0000 C   0  0  0  0  0  0
    0.8660   -1.2500    0.0000 C   0  0  0  0  0  0
    0.0000   -1.7500    0.0000 C   0  0  0  0  0  0
   -0.8660   -1.2500    0.0000 C   0  0  0  0  0  0
   -0.8660   -0.2500    0.0000 C   0  0  0  0  0  0
   -0.0000    1.2500    0.0000 C   0  0  0  0  0  0
    0.8660    1.7500    0.0000 O   0  0  0  0  0  0
   -0.8660    1.7500    0.0000 O   0  0  0  0  0  0
   -1.7321   -1.7500    0.0000 C   0  0  0  0  0  0
   -2.5981   -1.2500    0.0000 C   0  0  0  0  0  0
    1.7321   -1.7500    0.0000 C   0  0  0  0  0  0
    2.5981   -1.2500    0.0000 C   0  0  0  0  0  0
  1  21  0     0  0
  2  32  0     0  0
  3  41  0     0  0
  4  52  0     0  0
  5  61  0     0  0
  6  12  0     0  0
  1  71  0     0  0
  7  81  0     0  0
  7  92  0     0  0
  5 101  0     0  0
 10 111  1     0  0
  3 121  0     0  0
 12 131  6     0  0
M  END"""


@given("I have valid stereochem compounds within a ChemDraw file")
def step(context):
    # pull the contents of our chemdraw test file
    fn = os.path.join(os.path.dirname(__file__), 'files/behave_cdxml.cdxml')
    # convert the chemdraw file contents to mol
    mols = [mol.write("smi").split("\t")[0] for mol in readfile('cdxml', fn)]
    #file_contents = "".join(mols)
    print(len(mols))
    context.post_data["type"] = "Smiles"
    context.post_data["objects"] = mols
    # also populate our inchi list


@given("I have valid stereochem compounds within a SD file")
def step(context):
    # pull the contents of our SD test file
    fn = os.path.join(os.path.dirname(__file__), 'files/behave_sdf.sdf')
    mols = [mol.write("smi").split("\t")[0] for mol in readfile('sdf', fn)]
    print(len(mols))
    context.post_data["type"] = "Smiles"
    context.post_data["objects"] = mols


@then("I {action} my cbh_compound_batch to {projkey}")
def step(context, action=None, projkey=None, responsecode=202):
    from cbh_core_model.models import Project
    if action == "validate":
        # something
        path = "/dev/cbh_compound_batches/validate/"
        func = context.api_client.post
        context.post_data["projectKey"] = projkey
        resp = func(
            path,
            format='json',
            data=context.post_data,
        )
        assert int(resp.status_code) == int(responsecode)
    elif action == 'create':
        path = "/dev/cbh_compound_batches/"
        func = context.api_client.post
        context.post_data["projectKey"] = projkey
        # print(context.post_data)
        resp = func(
            path,
            format='json',
            data=context.post_data,
        )
        assert resp.status_code == 201


@then("I {action} my stereochem compounds to {projkey}")
def step(context, action=None, projkey=None, responsecode=202):
    from cbh_core_model.models import Project
    if action == "validate":
        # something
        path = "/dev/cbh_compound_batches/validate_list/"
        func = context.api_client.post
        context.post_data["projectKey"] = projkey
        resp = func(
            path,
            format='json',
            data=context.post_data,
        )
        context.post_data['current_batch'] = context.ser.deserialize(
            resp.content)["currentBatch"]
        assert int(resp.status_code) == int(responsecode)
    elif action == 'create':
        path = "/dev/cbh_compound_batches/multi_batch_save/"
        func = context.api_client.post
        context.post_data["projectKey"] = projkey
        # print(context.post_data)
        resp = func(
            path,
            format='json',
            data=context.post_data,
        )
        assert resp.status_code == 201


@then('retain its stereochemistry')
def step(context, action=None, projkey=None):
    from cbh_core_model.models import Project, CBHCompoundBatch
    from rdkit import Chem
    from rdkit.Chem import AllChem, inchi

    path = "/dev/cbh_compound_batches/"
    resp = context.api_client.get(
        path,
        format='json',
        data=context.post_data,
    )

    reg_cmpds = context.ser.deserialize(resp.content)["objects"]
    # retrieve registered inchi
    reg_inchi = reg_cmpds[0]['standardInchi']
    # convert our ctab mol to inchi
    #m = Chem.MolFromMolBlock(context.post_data["ctab"])
    #mol_inchi = inchi.MolToInchi(m)

    # we are now using a hard coded inchi from Chemicalize
    mol_inchi = context.inchi

    # assert they are equal
    assert mol_inchi == reg_inchi


@then("retain all their stereochemistry")
def step(context, action=None, projkey=None):
    # something here
    from cbh_core_model.models import Project, CBHCompoundBatch
    from rdkit import Chem
    from rdkit.Chem import AllChem, inchi

    path = "/dev/cbh_compound_batches/"
    resp = context.api_client.get(
        path,
        format='json',
        data=context.post_data,
    )

    reg_cmpds = context.ser.deserialize(resp.content)["objects"]
    reg_inchis = []
    # get a list of inchis from the reponse
    for cmpd in reg_cmpds:
        reg_inchis.append(cmpd['standardInchi'].strip())

    fn = os.path.join(os.path.dirname(__file__), 'files/inchi-list.txt')
    inchis = [mol.write("inchi").split("\t")[0].strip()
              for mol in readfile('inchi', fn)]

    # do an array subtraction of the hardcoded inchis from the registered inchis
    # print(set(inchis))
    print(len(inchis))
    # print(set(reg_inchis))
    print(len(reg_inchis))
    diff = list(set(inchis) - set(reg_inchis))
    print(len(diff))
    # print(diff)
    assert len(diff) == 0
