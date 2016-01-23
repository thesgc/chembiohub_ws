Feature: User can update projects with appropriate permissions


    Scenario: User cannot update projects if they have created them and logged out
        Given I create a project as before
        Given I take the first project in the list and change the name to Foo
        Given I log out
        When I patch the updated first project in the list back to the system
        Then the project update response is unauthorized



    Scenario: User cannot update projects if they have created them and removed all of their permissions
        Given I create a project as before
        Then I can list the projects on the system
        Given I remove all of the testusers permissions
        Given testuser has the cbh_core_model.add_project permission
        Given I take the first project in the list and change the name to Foo
        When I patch the updated first project in the list back to the system
        Then the project update response is unauthorized




    Scenario: User cannot update projects if they have created them and removed their ownership rights and have only project add rights
        Given I create a project as before
        Then I can list the projects on the system
        Given I remove all of the testusers permissions
        Given testuser has the cbh_core_model.add_project permission
        Given I take the first project in the list and change the name to Foo
        When I patch the updated first project in the list back to the system
        Then the project update response is unauthorized


    Scenario: User cannot update projects if they have created them and removed their ownership rights and have only editor rights
        Given I create a project as before
        Given I take the first project in the list and change the name to Foo
        Given I remove all of the testusers permissions
        Given I make testuser an editor of the first project in the list
        When I patch the updated first project in the list back to the system
        Then the project update response is unauthorized



    Scenario: User cannot update projects if they have created them and removed their ownership rights and have only viewer rights
        Given I create a project as before
        Given I remove all of the testusers permissions
        Given I make testuser a viewer of the first project in the list
        When I patch the updated first project in the list back to the system
        Then the project update response is unauthorized