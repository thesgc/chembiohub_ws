Feature: List and create projects

    Scenario: User can list the projects in ChemReg
        Given I have loaded the fixtures for project types
        Given testuser
        When I log in testuser
        Then I can list the projects on the system




