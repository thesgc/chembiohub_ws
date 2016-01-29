Feature: A user can save a compound with and without a structure and can clone a compound by remoiving the IDs and re submitting

    Scenario: A blinded compound batch or inventory record is created
        Given I have loaded the fixtures for project types and data types
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 4
        Given I create a project JSON by adding the saved search project type and the static custom field config schema for saved searches
        When I POST a project to cbh_projects
        Then the project is created
        When I refresh the user object
        Given I have a compound batch with no structure
        and I add the blinded batch id to my compound POST data as EMPTY_STRING
        and I add the project key to the compound data
        When I submit the compound by POST request
        and I reindex the compound
        Then a compound batch is created
        and the blinded batch ID is generated