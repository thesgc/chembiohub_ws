Feature: AssayReg data overview


    Scenario: A data form config can be linked to a project and then retrieved as a nested JSON
        Given I have created a data form config and a project as before and I list the projects
        and I take the first project in the list and I link it to my data form config via admin ui
        When I save the data form config via the admin ui
        Then All ancestor data form configs are present