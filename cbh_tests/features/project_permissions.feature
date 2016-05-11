Feature: User can get and update permissions on a project via the API if they are an owner

    
    Scenario: User can get permissions if owner
        Given I create a project as before
        Then I can list the projects on the system
        When I request the owner permissions for the first project in the list
        Then the response permission response is OK
        Then the permissions response shows that I am an owner

    
    Scenario: Owner cannot remove themselves
        Given I create a project as before
        When I request the owner permissions for the first project in the list 
        When I patch back the owner permissions as empty 
        Then the permissions response is unauthorized


    
    Scenario: Cannot add someone without the add project permission as owner
        Given I create a project as before
        Given anotheruser exists
        When I request the owner permissions for the first project in the list 
        When I patch back the owner permissions with anotheruser added
        Then the permissions response is unauthorized

    Scenario: Can add a superuser as owner as owner
        Given I create a project as before
        Given anotheruser exists
        Given I give anotheruser superuser permissions
        When I request the owner permissions for the first project in the list 
        When I patch back the owner permissions with anotheruser added
        Then the permissions response is accepted

    
    Scenario: Can add someone with the add project permission as owner
        Given I create a project as before
        Given anotheruser exists
        Given I give anotheruser add project permissions
        When I request the owner permissions for the first project in the list 
        When I patch back the owner permissions with anotheruser added
        Then the permissions response is accepted


    Scenario: Can add anyone as an editor regardless of if they have the add project permission
        Given I create a project as before
        Given anotheruser exists
        When I request the editor permissions for the first project in the list
        When I patch back the editor permissions with anotheruser added
        Then the permissions response is accepted