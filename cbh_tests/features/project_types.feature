Feature: Project Types provide project templates


    Scenario: User can list the project types
        Given I have loaded the fixtures for project types
        Given testuser
        When I log in testuser
        Then I can list the projects types on the system and there are 4


    Scenario: Default project type can be changed via the admin and there will always be only one
        Given I have loaded the fixtures for project types
        Given testuser
        When I log in testuser
        Then I can list the projects types on the system and there are 4
        Then the default project type has an id of 2
        When I take the project type model with id 1 and set it to default
        Then the previous default project type has been removed by the system and only the new one with id 1 remains


    Scenario: User can list the project types and there is one project type with the saved_search_project_type attribute and a static project template
        Given I have loaded the fixtures for project types
        Given testuser
        When I log in testuser
        Then I can list the projects types on the system and there are 4
        and the project type with id 3 has the saved_search_project_type
        and the project type with id 3 has a static project template to allow saved searches to be saved as compound batches

