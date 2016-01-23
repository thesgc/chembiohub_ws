Feature: Names of projects are unique and kept up to data with labels in django admin for permissions and 

    Scenario: Names of project can be changed
        Given I create a project as before
        Given I take the first project in the list and change the name to Foo
        When I patch the updated first project in the list back to the system
        Then project update response is accepted
        Then the name has changed to Foo and there is one project in the list

    Scenario: Names of the permissions, project key and custom field configs of a project are changed after name change
        Given I create a project as before
        Given I take the first project in the list and change the name to Foo
        When I patch the updated first project in the list back to the system
        Then project update response is accepted
        Then the project permission name matches the new name of the project
        Then the custom field config name matches the new name of the project
        Then the project key has changed


    Scenario: Project cannot be renamed to something that is already on the system
        Given I create a project as before
        Given I create a project JSON by adding one of these project types and some custom fields and Foo as a name
        When I POST a project to cbh_projects
        Then the project is created
        Then I can list the projects on the system
        Given I take the first project in the list and change the name to be equal to the second one
        When I patch the updated first project in the list back to the system
        Then the project update response is conflict


    Scenario: Project names must be unique when creating
        Given I create a project as before
        When I POST a project to cbh_projects
        Then The project is not created and there is a conflict

