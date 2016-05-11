Feature: I can search for molecules by substructure and flexmatch

    @wip
    Scenario: I start the qcluster
        Given I start the qcluster
        Given I register propane butane benzene and ethyl benzene via SMILES


        When I stop the qcluster
        Then I see the right qcluster output