Feature: List and create projects

    Scenario: User can list the projects in ChemReg
        Given testuser
        When I log in testuser
        Then I can list the projects on the system


    Scenario: User can list the project types on the system and these have been added by data migrations
        Given testuser
        When I log in testuser   
        Then I can list the projects types on the system and there are 3


