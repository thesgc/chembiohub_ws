Feature: List and create projects

    Scenario: User can list the projects in ChemReg
        Given testuser
        When I log in testuser
        Then I can list the projects on the system


    Scenario: User can list the project types on the system and these have been added by data migrations
        Given testuser
        When I log in testuser   
        Then I can list the projects types on the system and there are 3


    Scenario: User can create a project
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and a name
        When I POST a project to cbh_projects
        Then the project is created
        When I POST a project to cbh_projects



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