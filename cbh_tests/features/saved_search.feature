Feature: A user can save a search and owns it, organisations can also own a saved search


    Scenario: A saved search is created
        Given I have loaded the fixtures for project types
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 4
        Given I create a project JSON by adding the saved search project type and the static custom field config schema for saved searches
        When I POST a project to cbh_projects
        Then the project is created
        When I refresh the user object
        Given A URL to redirect the user to and a GET request URL and an alias and description for my saved search
        and I add the blinded batch id as EMPTY_ID
        When I send the search by POST request
        Then The saved search response is created


    Scenario: A user can list the compound batches on the system and does not see this saved search as it is in a saved search project
        Given I create a saved search as before
        When I list compound batches in the system
        Then the compound batch list response is OK
        and I see no compound batches


    Scenario: A user can list the saved searches on the system and sees this
        Given I create a saved search as before
        When I list saved searches in the system
        Then the saved search response is OK
        and I see my saved search




#Scenario also required for checking get detail on saved search API vs compound batch API but this is not implemented and not necessary for our tool