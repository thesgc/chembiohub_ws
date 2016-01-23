Feature: User can update projects with appropriate permissions

    Scenario: User can update projects if they have created them
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and a name
        When I POST a project to cbh_projects
        Then the project is created
        Then I can list the projects on the system
        Given I take the first project in the list and change the name
        When I patch the updated first project in the list back to the system
        Then project update response is accepted
        Then I can list the projects on the system
        Then the project name has changed




    Scenario: User cannot update projects if they have created them and logged out
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and a name
        When I POST a project to cbh_projects
        Then the project is created
        Then I can list the projects on the system
        Given I take the first project in the list and change the name
        Given I log out
        When I patch the updated first project in the list back to the system
        Then the project update response is unauthorized



    Scenario: User cannot update projects if they have created them and removed all of their permissions
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and a name
        When I POST a project to cbh_projects
        Then the project is created
        Then I can list the projects on the system
        Given I remove all of the testusers permissions
        Given testuser has the cbh_core_model.add_project permission
        Given I take the first project in the list and change the name
        When I patch the updated first project in the list back to the system
        Then the project update response is unauthorized




    Scenario: User cannot update projects if they have created them and removed their ownership rights and have only project add rights
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and a name
        When I POST a project to cbh_projects
        Then the project is created
        Then I can list the projects on the system
        Given I remove all of the testusers permissions
        Given testuser has the cbh_core_model.add_project permission
        Given I take the first project in the list and change the name
        When I patch the updated first project in the list back to the system
        Then the project update response is unauthorized


    Scenario: User cannot update projects if they have created them and removed their ownership rights and have only editor rights
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and a name
        When I POST a project to cbh_projects
        Then the project is created
        Then I can list the projects on the system
        Given I take the first project in the list and change the name
        Given I remove all of the testusers permissions
        Given I make testuser an editor of the first project in the list
        When I patch the updated first project in the list back to the system
        Then the project update response is unauthorized



    Scenario: User cannot update projects if they have created them and removed their ownership rights and have only viewer rights
        Given testuser
        Given testuser has the cbh_core_model.add_project permission
        When I log in testuser
        Then I can list the projects types on the system and there are 3
        Given I create a project JSON by adding one of these project types and some custom fields and a name
        When I POST a project to cbh_projects
        Then the project is created
        Then I can list the projects on the system
        Given I remove all of the testusers permissions
        Given I make testuser a viewer of the first project in the list
        When I patch the updated first project in the list back to the system
        Then the project update response is unauthorized