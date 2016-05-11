Feature: Restricted Fields

        Scenario: Data removed from output when field restricted
            Given testuser
            And I create a project as before
            And I add a restricted field
            Given I make testuser an owner of the first project in the list
            When I patch the updated first project in the list back to the system
            When a compound exists in the project with this restricted field
            Given I remove all of the testusers permissions
            Given I make testuser a viewer of the first project in the list
            When I list compound batches in the system
            Then the output data from search does not contain the restricted field


        Scenario: Field not present in schemata when restricted
            Given testuser
            And I create a project as before
            And I add a restricted field
            Given I make testuser an owner of the first project in the list
            When I patch the updated first project in the list back to the system
            Then I can list the projects on the system
            And the first project in the list contains the restricted field in the schemata
            Given I remove all of the testusers permissions
            Given I make testuser a viewer of the first project in the list
            Then I can list the projects on the system
            Then the first project in the list does not contain the restricted field in the schemata

