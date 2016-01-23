Feature: The required initial data can be loaded from fixtures or has been added by migration

    Scenario: User can list the project types on the system and these have been added by fixtures
        Given I have loaded the fixtures for project types and data types
        Given testuser
        When I log in testuser   
        Then I can list the projects types on the system and there are 3


    Scenario: User can list the data types on the system and these have been added by fixtures
        Given I have loaded the fixtures for project types and data types
        Given testuser
        When I log in testuser   
        Then I can list the data types on the system and there are 4
