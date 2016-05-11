Feature: The django qcluster is installed and can be started and stopped

    Scenario: I start the qcluster
        Given I start the qcluster
        When I stop the qcluster
        Then I see the right qcluster output