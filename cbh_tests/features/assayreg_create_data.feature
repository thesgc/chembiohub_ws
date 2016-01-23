Feature: AssayReg data overview


    Scenario: A data form config can be linked to a project and then retrieved as a nested JSON
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and Bar as a name
        When I POST a project to cbh_projects
        Then the project is created
        Given I create new custom field configs and data form configs based on the data given
        When I save the data form config via the admin ui
        Then All ancestor data form configs are present 
        Then I can list the projects on the system
        Given I take the first project in the list and I link it to my data form config via admin ui
        When I get my single project via the data form config API as in assayreg data overview
        Then there is a nest of data form configs down to l3

