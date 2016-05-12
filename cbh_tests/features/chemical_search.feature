Feature: I can search for molecules by substructure and flexmatch

    @wip
    Scenario: I start the qcluster
        Given I start the qcluster


        When I stop the qcluster
        Then I see the right qcluster output