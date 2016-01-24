Feature: AssayReg data overview


    Scenario: A data form config can be linked to a project and then retrieved as a nested JSON
        Given I have created a data form config and a project as before and I list the projects
        and I take the first project in the list and I link it to my data form config via admin ui
        When I save the data form config via the admin ui
        Then All ancestor data form configs are present



    Scenario: List projects with forms
        Given I set up a project and data form config as before
        When I list the projects with forms
        Then the first project in the list has the expected data form config



    Scenario: Data overview page first save of project information
        Given I make a request to projects with forms as before
        When I add data to the l0 data point classification template 
        And I POST it to the data point classification resource
        Then the l0 data point classification is created for my project


    Scenario: Data overview page can list first data point classification
        Given I create a datapoint classification as before
        When I list the nested datapoint classifications in the project
        Then the nested classification response is as expected and the resource URI is ready
