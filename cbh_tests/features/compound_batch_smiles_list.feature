Feature: A user can register a list of SMILES strings

    Scenario: A user can validate a SMILES list
        Given I start the qcluster
        Given I set up the SMILES data
        When I validate propane butane benzene and ethyl benzene via SMILES
        Then the response from post validate list is accepted
        When I stop the qcluster
        Then I see the right qcluster output


    Scenario: A user can save a SMILES list
        Given I start the qcluster
        Given I set up the SMILES data
        When I validate propane butane benzene and ethyl benzene via SMILES
        Then the response from post validate list is accepted
        When I take the response from post validate drawn and post it to multi batch save
        then the response from multi batch save is created
        When I list compound batches in the system
        Then 4 compound batches have been created
        When I stop the qcluster
        Then I see the right qcluster output

