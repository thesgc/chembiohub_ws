Feature: I can create datapointclassifications


    Scenario: List the data form configs for a project
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and a name
        When I POST a project to cbh_projects
        Then the project is created
        Given I create new custom field configs and data form configs based on the data given
        Given I link my project to my data form config


