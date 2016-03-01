Feature: Project List


    Scenario Outline: Project list not logged in
        Given a User
        and my user is member of a group
        and I have a valid molfile
        and a valid project exists proja
        and I automatically have editor permissions as creator
        and I remove my permissions
        and I have<editor> given <me_or_group> editor privileges for proja
        and I have<viewer> given <me_or_group> viewer privileges for proja        
        When do not log in
        When I get projects list the response code will be 401



    Scenario Outline: Project list includes my project if I have privileges
        Given a User
        and my user is member of a group
        and I have a valid molfile
        and a valid project exists proja
        and I automatically have editor permissions as creator
        and I remove my permissions
        and I have<editor> given <me_or_group> editor privileges for proja
        and I have<viewer> given <me_or_group> viewer privileges for proja        
        When <dologin> log in
        Then I <action> projects and the <responsecode>


        Examples: Validation api
        |me_or_group           |   dologin    |   editor  |   viewer  |   responsecode   | action |
        |myself|      I         |           |           |  response code will be 200 with proja in it   | get |
        |myself|       I        |           |   nt     |  response code will be 200 with proja in it    | get |
        |myself|        I       |   nt     |           |  response code will be 200 with proja in it    | get |
        |myself|   I do not      |           |           |  response code will be 401            | get |
                |myself|         I      |   nt     |   nt     |  response code will be 200 without proja in it  | get |
