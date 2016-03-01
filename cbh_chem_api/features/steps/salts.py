from behave import given, when, then
import json
import os
from django.contrib.auth.models import User, Group
from pybel import readfile


@given('I have valid salted compound')
def step(context):
    context.inchi = "InChI=1S/2CHF3O3S.Cu/c2*2-1(3,4)8(5,6)7;/h2*(H,5,6,7);/q;;+2/p-2"
    context.post_data["ctab"] = """


 17 14  0  0000  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Cu  0  2  0  0  0  0  0  0  0  0  0
    2.8286   -1.2375    0.0000 O   0  5  0  0  0  0  0  0  0  0  0
    2.8286   -0.4125    0.0000 S   0  0  0  0  0  0  0  0  0  0  0
    3.6536   -0.4125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0
    2.0036   -0.4125    0.0000 O   0  0  0  0  0  0  0  0  0  0  0
    2.8286    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    2.8286    1.2375    0.0000 F   0  0  0  0  0  0  0  0  0  0  0
    3.6536    0.4125    0.0000 F   0  0  0  0  0  0  0  0  0  0  0
    2.0036    0.4125    0.0000 F   0  0  0  0  0  0  0  0  0  0  0
    0.0000   -5.3625    0.0000 O   0  5  0  0  0  0  0  0  0  0  0
    0.0000   -4.5375    0.0000 S   0  0  0  0  0  0  0  0  0  0  0
    0.8250   -4.5375    0.0000 O   0  0  0  0  0  0  0  0  0  0  0
   -0.8250   -4.5375    0.0000 O   0  0  0  0  0  0  0  0  0  0  0
    0.0000   -3.7125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0
    0.0000   -2.8875    0.0000 F   0  0  0  0  0  0  0  0  0  0  0
    0.8250   -3.7125    0.0000 F   0  0  0  0  0  0  0  0  0  0  0
   -0.8250   -3.7125    0.0000 F   0  0  0  0  0  0  0  0  0  0  0
  2  3  1  0
  3  4  2  0
  3  5  2  0
  3  6  1  0
  6  7  1  0
  6  8  1  0
  6  9  1  0
 10 11  1  0
 11 12  2  0
 11 13  2  0
 11 14  1  0
 14 15  1  0
 14 16  1  0
 14 17  1  0
M  CHG  3   1   2   2  -1  10  -1
M  END"""


@given("I have valid salted compounds within a ChemDraw {format} file")
def step(context, format=None):
    # pull the contents of our chemdraw test file
    fp = 'files/behave-triflate.%s' % (format)
    fn = os.path.join(os.path.dirname(__file__), fp)
    # convert the chemdraw file contents to mol
    print(fp)
    print(fn)
    print(format)
    mols = [mol.write("smi").split("\t")[0] for mol in readfile(format, fn)]
    #file_contents = "".join(mols)
    print(len(mols))
    context.post_data["type"] = "Smiles"
    context.post_data["objects"] = mols
    # also populate our inchi list


@then('retain its salt')
def step(context, action=None, projkey=None):
    from cbh_chembl_model_extension.models import CBHCompoundBatch
    from cbh_core_model.models import Project

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
