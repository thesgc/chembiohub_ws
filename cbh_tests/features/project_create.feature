Feature: user can create projects


    Scenario: User can create a project if they have permission to
        Given I have loaded the fixtures for project types and data types
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 4
        Given I create a project JSON by adding one of these project types and some custom fields and Bar as a name
        When I POST a project to cbh_projects
        Then the project is created




    Scenario: User cannot create a project if they do not have that permission
        Given I have loaded the fixtures for project types and data types
        Given testuser
        When I log in testuser
        Then I can list the projects types on the system and there are 4
        Given I create a project JSON by adding one of these project types and some custom fields and Bar as a name
        When I POST a project to cbh_projects
        Then the project is created not as unauthorized


    Scenario: User becomes owner of a project that they create
        Given I create a project as before
        Then the project creator is automatically an owner
