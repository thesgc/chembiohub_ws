Feature: user can create projects


    Scenario: User can create a project if they have that permission
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and a name
        When I POST a project to cbh_projects
        Then the project is created


    Scenario: User cannot create a project if they have that permission
        Given testuser
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and a name
        When I POST a project to cbh_projects
        Then the project is created not as unauthorized


    Scenario: Project names must be unique
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and a name
        When I POST a project to cbh_projects
        Then the project is created
        When I POST a project to cbh_projects
        Then The project is not created and there is a conflict


