Feature: I can search for molecules by substructure and flexmatch

    Scenario: I can search for an exact match with structure search
        Given I start the qcluster
        Given I save propane butane benzene and ethyl benzene via SMILES
        When I POST the molfile for propane as a flexmatch
        Then I get a chemical search id in the response
        When I run a chemical search
        Then I get 1 chemical search result
        When I stop the qcluster
        Then I see the right qcluster output


    Scenario: I can search for a substructure match with structure search
        Given I start the qcluster
        Given I save propane butane benzene and ethyl benzene via SMILES
        When I POST the molfile for propane as a with_substructure
        Then I get a chemical search id in the response
        When I run a chemical search
        Then I get 3 chemical search results
        When I stop the qcluster
        Then I see the right qcluster output


