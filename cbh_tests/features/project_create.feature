Feature: user can create projects


    Scenario: User can create a project if they have permission to
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and Bar as a name
        When I POST a project to cbh_projects
        Then the project is created




    Scenario: User cannot create a project if they do not have that permission
        Given testuser
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and Bar as a name
        When I POST a project to cbh_projects
        Then the project is created not as unauthorized


    Scenario: User becomes owner of a project that they create
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and Bar as a name
        When I POST a project to cbh_projects
        Then I can list the projects on the system
        Then the project is created
        Then the project creator is automatically an owner




    Scenario: User can create a project with a data form config
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and Bar as a name
        Given I also add 2 data form configs to that project
        When I POST a project to cbh_projects_with_forms
        Then I can list the projects on the system
        Then the project is created