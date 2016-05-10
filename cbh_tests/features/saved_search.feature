Feature: A user can save a search and owns it, organisations can also own a saved search

    @wipnot
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
        and I add the project key
        When I send the search by POST request
        Then The saved search response is created
        When I reindex the saved search 
        Then The saved search index response is OK


    @wipnot
    Scenario: A user can list the compound batches on the system and does not see this saved search as it is in a saved search project
        Given I create a saved search as before
        When I list compound batches in the system with get_list_elasticsearch
        Then the get list elasticsearch response is OK
        and I see no compound batches


    @wipnot
    Scenario: A user can list the saved searches on the system and sees this
        Given I create a saved search as before
        When I list saved searches in the system with get_list_elasticsearch
        Then the get list elasticsearch response is OK
        and I see my saved search